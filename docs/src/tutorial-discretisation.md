# [Discretisation methods](@id tutorial-discretisation-methods)

In this tutorial, we will explore the different discretisation methods available in the OptimalControl.jl package.
These methods are used to convert a continuous-time optimal control problem (OCP) into a discrete-time nonlinear programming (NLP) problem, which can then be solved using numerical optimization techniques.

Let us import the necessary packages and define the optimal control problem ([Goddard problem](@ref tutorial-goddard)) we will use as an example throughout this tutorial.

```@example main-disc
using BenchmarkTools  # for benchmarking
using DataFrames      # to store the results
using OptimalControl  # to define the optimal control problem and more
using NLPModels       # to retrieve information about the NLP problem 
using NLPModelsIpopt  # to solve the problem via a direct method
using Plots           # to plot the solution

t0 = 0      # initial time
r0 = 1      # initial altitude
v0 = 0      # initial speed
m0 = 1      # initial mass
vmax = 0.1  # maximal authorized speed
mf = 0.6    # final mass to target

# Goddard problem function
function goddard_problem()
    ocp = @def begin # definition of the optimal control problem

    tf ∈ R, variable
    t ∈ [t0, tf], time
    x = (r, v, m) ∈ R³, state
    u ∈ R, control

    x(t0) == [r0, v0, m0]
    m(tf) == mf
    0 ≤ u(t) ≤ 1
    r(t) ≥ r0
    0 ≤ v(t) ≤ vmax

    ẋ(t) == F0(x(t)) + u(t) * F1(x(t))

    r(tf) → max

    end
    return ocp
end

ocp = goddard_problem()

# Dynamics
Cd = 310
Tmax = 3.5
β = 500
b = 2

F0(x) = begin
    # Uncontrolled dynamics: gravity and drag
    r, v, m = x
    D = Cd * v^2 * exp(-β*(r - 1)) # Aerodynamic drag
    return [ v, -D/m - 1/r^2, 0 ]   # [dr/dt, dv/dt, dm/dt]
end

F1(x) = begin
    # Control dynamics: thrust contribution
    r, v, m = x
    return [ 0, Tmax/m, -b*Tmax ]   # [dr/dt, dv/dt, dm/dt] due to thrust
end
nothing # hide
```

## Discretisation schemes

When calling `solve`, the option `scheme=...` can be used to specify the discretization scheme. In addition to the default implicit `:midpoint` method, other available options include the implicit `:trapeze` method (also known as Crank–Nicolson), and Gauss–Legendre collocation schemes with 2 and 3 stages: `:gauss_legendre_2` and `:gauss_legendre_3`, of order 4 and 6 respectively. Note that higher-order methods typically result in larger NLP problems for the same number of time steps, and their accuracy also depends on the smoothness of the problem.

Let us first solve the problem with the default `:midpoint` method and display the solution.

```@example main-disc
sol = solve(ocp; scheme=:midpoint, display=false)
@assert successful(sol) "Solution failed"
plot(sol; size=(800, 800))
```

Let us now compare different discretization schemes to evaluate their accuracy and performance. See the [Strategy options](@id manual-solve-strategy-options) for information on default values for parameters like `tol`.

```@example main-disc
# Solve the problem with different discretization methods
solutions = []
data = DataFrame(
    Scheme=Symbol[],
    Time=Float64[],
    Objective=Float64[], 
    Iterations=Int[]
)
schemes = [
    :euler, 
    :euler_implicit,
    :trapeze, 
    :midpoint, 
    :gauss_legendre_2, 
    :gauss_legendre_3
]
for scheme in schemes
    bt = @btimed solve($ocp; scheme=$scheme, tol=1e-8, display=false)
    local sol = bt.value
    #@assert successful(sol) "Solution failed for scheme=$scheme"
    push!(solutions, (scheme, sol))
    push!(data, (
        Scheme=scheme, 
        Time=bt.time, 
        Objective=objective(sol), 
        Iterations=iterations(sol)
    ))
end
println(data)
```

```@example main-disc
# Plot the results
x_style = (legend=:none,)
p_style = (legend=:none,)
u_style = (legend=:topright,)
styles = (state_style=x_style, costate_style=p_style, control_style=u_style)

plt = plot()
for (i, (scheme, sol)) in enumerate(solutions)
    plt = plot!(sol; label=string(scheme), color=i, styles...)
end
plot(plt; size=(800, 800))
```

## Large problems and AD backend

For some large problems, you may notice that the solver takes a long time before the iterations actually begin. This is due to the computation of sparse derivatives — specifically, the Jacobian of the constraints and the Hessian of the Lagrangian — which can be quite costly. One possible alternative is to set the option `backend=:manual` in the modeler, which uses simpler sparsity patterns. The resulting matrices are faster to compute but are also less sparse, so this represents a trade-off between automatic differentiation (AD) preparation time and the efficiency of the optimization itself.

```@example main-disc
disc = OptimalControl.Collocation(grid_size=1000, scheme=:gauss_legendre_3)
mod = OptimalControl.ADNLP(backend=:manual)
sol = OptimalControl.Ipopt()
solve(ocp; discretizer=disc, modeler=mod, solver=sol)
nothing # hide
```

This explicit approach with strategy instances is less high-level than the descriptive method used earlier. It provides more control over the solving process and is detailed in the [explicit mode documentation](@extref OptimalControl manual-solve-explicit).

Let us now compare the performance of the two backends. We will use the `@btimed` macro from the `BenchmarkTools` package to measure the time taken for both the preparation of the NLP problem and the execution of the solver. Separating these measurements allows us to understand whether the bottleneck is in the discretization/modeling phase (which depends on the AD backend) or in the solver phase (which depends on the sparsity of the matrices). We use a lower-level approach with explicit strategy instances to separately measure the discretization/modeling time and the solver time.

```@example main-disc
# DataFrame to store the results
data = DataFrame(
    Backend=Symbol[],
    NNZO=Int[],
    NNZJ=Int[],
    NNZH=Int[],
    PrepTime=Float64[],
    ExecTime=Float64[],
    Objective=Float64[], 
    Iterations=Int[]
)

# The different AD backends to test
backends = [:optimized, :manual]

# Loop over the backends
for adnlp_backend in backends

    # Create strategy instances
    discretizer = OptimalControl.Collocation(grid_size=1000, scheme=:gauss_legendre_3)
    modeler = OptimalControl.ADNLP(backend=adnlp_backend)
    solver = OptimalControl.Ipopt(print_level=0)

    # Initial guess
    init = build_initial_guess(ocp, nothing)

    # Discretize and build NLP model
    # This lower-level approach with discretize/nlp_model/solve is detailed in the [tutorial on direct transcription](@ref tutorial-nlp)
    bt = @btimed begin
        docp = discretize($ocp, $discretizer)
        nlp_model(docp, $init, $modeler)
    end
    prepa_time = bt.time
    nlp = bt.value

    # Get the number of non-zero elements
    nnzo = get_nnzo(nlp) # Gradient of the Objective
    nnzj = get_nnzj(nlp) # Jacobian of the constraints
    nnzh = get_nnzh(nlp) # Hessian of the Lagrangian

    # Solve the problem
    bt = @btimed solve($nlp, $solver; display=false)
    exec_time = bt.time
    nlp_sol = bt.value

    # Store the results in the DataFrame
    push!(data, 
        (
            Backend=adnlp_backend, 
            NNZO=nnzo, 
            NNZJ=nnzj, 
            NNZH=nnzh, 
            PrepTime=prepa_time, 
            ExecTime=exec_time, 
            Objective=nlp_sol.objective, 
            Iterations=nlp_sol.iter
        )
    )
end
println(data)
```

## Explicit time grid

The option `time_grid=...` allows you to provide the full time grid vector `t0, t1, ..., tf`, which is especially useful if a non-uniform grid is desired. In the case of a free initial and/or final time, you should provide a normalized grid ranging from 0 to 1. Note that `time_grid` overrides `grid_size` if both options are specified.

```@example main-disc
sol = solve(ocp; time_grid=[0, 0.1, 0.5, 0.9, 1], display=false)
println(time_grid(sol))
```
