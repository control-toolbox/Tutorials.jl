using OptimalControl  # to define the optimal control problem and more
using NLPModelsIpopt  # to solve the problem via a direct method
using Plots           # to plot the solution

t0 = 0      # initial time
r0 = 1      # initial altitude
v0 = 0      # initial speed
m0 = 1      # initial mass
vmax = 0.1  # maximal authorized speed
mf = 0.6    # final mass to target

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
Cd = 310
Tmax = 3.5
β = 500
b = 2

F0(x) = begin
    r, v, m = x
    D = Cd * v^2 * exp(-β*(r - 1)) # Drag force
    return [v, -D/m - 1/r^2, 0]
end

F1(x) = begin
    r, v, m = x
    return [0, Tmax/m, -b*Tmax]
end

# Solve the problem with different discretization methods
solutions = []
schemes = [
    :trapeze, :midpoint, :euler, :euler_implicit, :gauss_legendre_2, :gauss_legendre_3
]
for scheme in schemes
    sol = solve(ocp; disc_method=scheme, tol=1e-8, display=false)
    push!(solutions, (scheme, sol))
end

# Plot the results
x_style = (legend=:none,)
p_style = (legend=:none,)
styles = (state_style=x_style, costate_style=p_style)

scheme, sol = solutions[1]
plt = plot(sol; label=string(scheme), styles...)
for (scheme, sol) in solutions[2:end]
    plt = plot!(sol; label=string(scheme), styles...)
end
plot(plt; size=(800, 800))

# Benchmark the problem with different AD backends
using NLPModels
using DataFrames
using BenchmarkTools

# DataFrame to store the results
data = DataFrame(;
    Backend=Symbol[],
    NNZO=Int[],
    NNZJ=Int[],
    NNZH=Int[],
    PrepTime=Float64[],
    ExecTime=Float64[],
    Objective=Float64[],
    Iterations=Int[],
)

# The different AD backends to test
backends = [:optimized, :manual]

# Loop over the backends
for adnlp_backend in backends
    println("Testing backend: ", adnlp_backend)

    # Discretize the problem with a large grid size and Gauss-Legendre method
    bt = @btimed direct_transcription(
        $ocp; disc_method=:euler, grid_size=1000, adnlp_backend=($adnlp_backend)
    )
    docp, nlp = bt.value
    prepa_time = bt.time

    # Get the number of non-zero elements
    nnzo = get_nnzo(nlp) # Gradient of the Objective
    nnzj = get_nnzj(nlp) # Jacobian of the constraints
    nnzh = get_nnzh(nlp) # Hessian of the Lagrangian

    # Solve the problem
    bt = @btimed ipopt($nlp; print_level=0, mu_strategy="adaptive", tol=1e-8, sb="yes")
    nlp_sol = bt.value
    exec_time = bt.time

    # Store the results in the DataFrame
    push!(
        data,
        (
            Backend=adnlp_backend,
            NNZO=nnzo,
            NNZJ=nnzj,
            NNZH=nnzh,
            PrepTime=prepa_time,
            ExecTime=exec_time,
            Objective=(-nlp_sol.objective),
            Iterations=nlp_sol.iter,
        ),
    )
end
println(data)
