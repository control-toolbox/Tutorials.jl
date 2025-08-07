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

const t0 = 0      # initial time
const r0 = 1      # initial altitude
const v0 = 0      # initial speed
const m0 = 1      # initial mass
const vmax = 0.1  # maximal authorized speed
const mf = 0.6    # final mass to target

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

end;

# Dynamics
const Cd = 310
const Tmax = 3.5
const β = 500
const b = 2

F0(x) = begin
    r, v, m = x
    D = Cd * v^2 * exp(-β*(r - 1)) # Drag force
    return [ v, -D/m - 1/r^2, 0 ]
end

F1(x) = begin
    r, v, m = x
    return [ 0, Tmax/m, -b*Tmax ]
end
nothing # hide
```

## Discretisation schemes

When calling `solve`, the option `disc_method=...` can be used to specify the discretization scheme. In addition to the default implicit `:trapeze` method (also known as Crank–Nicolson), other available options include the implicit `:midpoint` method, and Gauss–Legendre collocation schemes with 2 and 3 stages: `:gauss_legendre_2` and `:gauss_legendre_3`, of order 4 and 6 respectively. Note that higher-order methods typically result in larger NLP problems for the same number of time steps, and their accuracy also depends on the smoothness of the problem.

Let us first solve the problem with the default `:trapeze` method and display the solution.

```@example main-disc
sol = solve(ocp; disc_method=:trapeze, display=false)
plot(sol; size=(800, 800))
```

Let us now compare different discretization schemes to evaluate their accuracy and performance.

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
    bt = @btimed solve($ocp; disc_method=$scheme, tol=1e-8, display=false)
    local sol = bt.value
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

scheme, sol = solutions[1]
plt = plot(sol; label=string(scheme), styles...)
for (scheme, sol) in solutions[2:end]
    plt = plot!(sol; label=string(scheme), styles...)
end
plot(plt; size=(800, 800))
```

## Large problems and AD backend

For some large problems, you may notice that the solver takes a long time before the iterations actually begin. This is due to the computation of sparse derivatives — specifically, the Jacobian of the constraints and the Hessian of the Lagrangian — which can be quite costly. One possible alternative is to set the option `adnlp_backend=:manual`, which uses simpler sparsity patterns. The resulting matrices are faster to compute but are also less sparse, so this represents a trade-off between automatic differentiation (AD) preparation time and the efficiency of the optimization itself.

```@example main-disc
solve(ocp; disc_method=:gauss_legendre_3, grid_size=1000, adnlp_backend=:manual)
nothing # hide
```

Let us now compare the performance of the two backends. We will use the `@btimed` macro from the `BenchmarkTools` package to measure the time taken for both the preparation of the NLP problem and the execution of the solver. Besides, we will also collect the number of non-zero elements in the Jacobian and Hessian matrices, which can be useful to understand the sparsity of the problem, thanks to the functions `get_nnzo`, `get_nnzj`, and `get_nnzj` from the `NLPModels` package. The problem is first discretized with the `direct_transcription` method and then solved with the `ipopt` solver, see the [tutorial on direct transcription](@ref tutorial-nlp) for more details.

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

    # Discretize the problem with a large grid size and Gauss-Legendre method
    bt = @btimed direct_transcription($ocp; 
        disc_method=:gauss_legendre_3, 
        grid_size=1000, 
        adnlp_backend=$adnlp_backend,
    )
    docp = bt.value
    nlp = model(docp)
    prepa_time = bt.time

    # Get the number of non-zero elements
    nnzo = get_nnzo(nlp) # Gradient of the Objective
    nnzj = get_nnzj(nlp) # Jacobian of the constraints
    nnzh = get_nnzh(nlp) # Hessian of the Lagrangian

    # Solve the problem
    bt = @btimed ipopt($nlp; 
        print_level=0, 
        mu_strategy="adaptive", 
        tol=1e-8,
        sb="yes",
    )
    nlp_sol = bt.value
    exec_time = bt.time

    # Store the results in the DataFrame
    push!(data, 
        (
            Backend=adnlp_backend, 
            NNZO=nnzo, 
            NNZJ=nnzj, 
            NNZH=nnzh, 
            PrepTime=prepa_time, 
            ExecTime=exec_time, 
            Objective=-nlp_sol.objective, 
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
