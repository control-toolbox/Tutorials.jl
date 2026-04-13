# Discrete continuation

By using the warm start option, it is easy to implement a basic discrete continuation method, in which a sequence of problems is solved by using each solution as the initial guess for the next problem. This approach typically leads to faster and more reliable convergence than solving each problem with the same initial guess and is particularly useful for problems that require a good initial guess to converge.

## Continuation on parametric OCP

The most concise way to perform discrete continuation is to define a function that returns the optimal control problem for a given value of the continuation parameter, and then solve a sequence of such problems.
We illustrate this using a simple double integrator problem, where the fixed final time is gradually increased.

First we load the required packages:

```@example main-cont
using DataFrames
using OptimalControl
using NLPModelsIpopt
using Printf
using Plots
```

The `init` parameter of the `solve` function allows providing an initial guess. See the [initial guess documentation](@extref OptimalControl manual-initial-guess) for more details.

We write a function that returns the OCP for a given final time:

```@example main-cont
function problem(T)

    ocp = @def begin

        t ∈ [0, T], time
        x ∈ R², state
        u ∈ R, control

        q = x₁
        v = x₂

        q(0) == 0
        v(0) == 0
        q(T) == 1
        v(T) == 0
        ẋ(t) == [v(t), u(t)]

        ∫(u(t)^2) → min

    end

    return ocp
end
nothing # hide
```

Then we perform the continuation with a simple *for* loop, using each solution to initialize the next problem. We wrap the continuation in a function to avoid global variables.

```@example main-cont
function continuation_parametric(T_range; init=nothing, scheme=:midpoint, grid_size=200)
    data = DataFrame(T=Float64[], Objective=Float64[], Iterations=Int[])
    for T ∈ T_range
        ocp = problem(T)
        sol = solve(ocp; init=init, display=false, scheme=scheme, grid_size=grid_size)
        @assert successful(sol) "Solution failed for T=$T"
        init = sol
        push!(data, (T=T, Objective=objective(sol), Iterations=iterations(sol)))
    end
    return data
end

data = continuation_parametric(range(1, 2, length=5))
println(data)
```

We can visualize the evolution of the objective and the number of iterations with respect to the final time.

```@example main-cont
plt_obj = plot(data.T, data.Objective;
    seriestype=:scatter,
    title="Double integrator",
    label="Objective",
    xlabel="Final time T",
    ylabel="∫u(t)² dt")

plt_iter = plot(data.T, data.Iterations;
    seriestype=:scatter,
    label="Iterations",
    xlabel="Final time T",
    ylabel="Number of iterations")

layout = grid(2, 1, heights=[0.5, 0.5])
plot(plt_obj, plt_iter; layout=layout, size=(800, 600))
```

## Continuation on parameter

As a second example, we solve a Goddard problem with a decreasing maximum thrust. We define a function that returns the optimal control problem for a given value of `Tmax`, similar to the first example.

Let us first define the Goddard problem. Note that the formulation below illustrates all types of constraints, and the problem could be written more compactly.

```@example main-cont
# Parameters
r0 = 1
v0 = 0
m0 = 1
mf = 0.6
x0 = [r0, v0, m0]
vmax = 0.1

# Dynamics
function F0(x)
    # Uncontrolled dynamics: gravity and drag
    r, v, m = x
    D = Cd * v^2 * exp(-β*(r - 1))  # Aerodynamic drag
    return [ v, -D/m - 1/r^2, 0 ]   # [dr/dt, dv/dt, dm/dt]
end
function F1(x, Tmax)
    # Control dynamics: thrust contribution
    r, v, m = x
    return [ 0, Tmax/m, -b*Tmax ]   # [dr/dt, dv/dt, dm/dt] due to thrust
end

# Parameters for the dynamics
Cd = 310
β = 500
b = 2
Tmax_0 = 3
Tmax_f = 1

# Goddard problem function that takes Tmax as parameter
function goddard_problem(Tmax)
    ocp = @def begin

        tf ∈ R, variable
        t ∈ [0, tf], time
        x ∈ R^3, state
        u ∈ R, control

        0.01 ≤ tf ≤ Inf

        r = x[1]
        v = x[2]
        m = x[3]
        x(0) == x0
        m(tf) == mf
        r0 ≤ r(t) ≤ r0 + 0.1
        v0 ≤ v(t) ≤ vmax
        mf ≤ m(t) ≤ m0
        0 ≤ u(t) ≤ 1
        ẋ(t) == F0(x(t)) + u(t) * F1(x(t), Tmax)

        r(tf) → max

    end
    return ocp
end

# Solve the problem with a reference value of Tmax
sol0 = solve(goddard_problem(Tmax_0); display=false, scheme=:midpoint, grid_size=200)
@printf("Objective for reference solution: %.6f\n", objective(sol0))
```

Then, we perform the continuation on the maximal thrust. We wrap the continuation in a function that redefines the OCP at each step with the new `Tmax` value.

```@example main-cont
function continuation_goddard(Tmax_range; init=nothing, scheme=:midpoint, grid_size=200)
    data = DataFrame(Tmax=Float64[], Objective=Float64[], Iterations=Int[])
    sols = Vector{Any}()
    sol = init
    for Tmax_local ∈ Tmax_range
        ocp = goddard_problem(Tmax_local)
        sol = solve(ocp; init=sol, display=false, scheme=scheme, grid_size=grid_size)
        @assert successful(sol) "Solution failed for Tmax=$Tmax_local"
        push!(data, (Tmax=Tmax_local, Objective=objective(sol), Iterations=iterations(sol)))
        push!(sols, sol)
    end
    return data, sols
end

data, sols = continuation_goddard(range(Tmax_0, Tmax_f, length=9); init=sol0)
println(data)
```

We plot now the objective with respect to the maximal thrust, as well as the solutions for `Tmax=3`, `Tmax=2`, and `Tmax=1`. The time is normalized for the solution plots to compare trajectories with different final times.

```@example main-cont
using Plots.PlotMeasures # for leftmargin

plt_obj = plot(data.Tmax, data.Objective;
    seriestype=:scatter,
    title="Goddard problem",
    label="r(tf)",
    xlabel="Maximal thrust (Tmax)",
    ylabel="Maximal altitude r(tf)")

plt_iter = plot(data.Tmax, data.Iterations;
    seriestype=:scatter,
    label="Iterations",
    xlabel="Maximal thrust (Tmax)",
    ylabel="Number of iterations")

# Find indices closest to desired Tmax values
Tmax_values = [1.0, 2.0, 3.0]
indices = [argmin(abs.(data.Tmax .- T)) for T in Tmax_values]

layout_metrics = grid(2, 1, heights=[0.5, 0.5])
plot(plt_obj, plt_iter; layout=layout_metrics, size=(800, 600), leftmargin=5mm)
```

We now plot the solutions for the selected `Tmax` values, using normalized time to compare trajectories with different final times.

```@example main-cont
plt_sol = plot()
for (i, idx) in enumerate(indices)
    plot!(plt_sol, sols[idx]; label="Tmax=$(round(data.Tmax[idx], digits=2))", time=:normalize, color=i)
end
plot(plt_sol; size=(800, 800), leftmargin=5mm)
```

## Best practices and limitations

The examples above illustrate a basic discrete continuation method. Here are some important considerations:

- **Step size**: The continuation step size should be chosen carefully. Too large a step may lead to convergence issues, while too small a step may be inefficient. Adaptive step size strategies can improve robustness but are not covered in this tutorial.

- **Convergence monitoring**: Always check that each solution converges successfully using `successful(sol)`. If a solution fails, consider reducing the step size or providing a better initial guess.

- **Advanced methods**: For more challenging problems, homotopic continuation methods with path following algorithms can be used. These methods adapt the continuation parameter automatically and can handle bifurcations. Such advanced techniques are beyond the scope of this tutorial.

- **Alternative approaches**: In some cases, it may be more efficient to directly parameterize the continuation parameter within the optimal control problem rather than solving a sequence of problems. OptimalControl plans to add this feature in a future release.
