# Discrete continuation

By using the warm start option, it is easy to implement a basic discrete continuation method, in which a sequence of problems is solved by using each solution as the initial guess for the next problem.
This approach typically leads to faster and more reliable convergence than solving each problem with the same initial guess and is particularly useful for problems that require a good initial guess to converge.

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

and write a function that returns the OCP for a given final time:

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

Then we perform the continuation with a simple *for* loop, using each solution to initialize the next problem.

```@example main-cont
init = ()
data = DataFrame(T=Float64[], Objective=Float64[], Iterations=Int[])
for T ∈ range(1, 2, length=5)
    ocp = problem(T) 
    sol = solve(ocp; init=init, display=false)
    global init = sol
    push!(data, (T=T, Objective=objective(sol), Iterations=iterations(sol)))
end
println(data)
```

## Continuation on global variable

As a second example, we show how to avoid redefining a new optimal control problem at each step by modifying the original one instead. More precisely, we solve a Goddard problem with a decreasing maximum thrust. By storing the value of `Tmax` in a global variable, we can simply update this variable and reuse the same problem throughout the continuation.

Let us first define the Goddard problem. Note that the formulation below illustrates all types of constraints, and the problem could be written more compactly.

```@example main-cont
# Parameters
r0 = 1
v0 = 0
m0 = 1
mf = 0.6
x0 = [r0, v0, m0]
vmax = 0.1

# Goddard problem definition
@def goddard begin

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
    ẋ(t) == F0(x(t)) + u(t) * F1(x(t))

    r(tf) → max

end

# Dynamics
function F0(x)
    r, v, m = x
    D = Cd * v^2 * exp(-β*(r - 1))
    return [ v, -D/m - 1/r^2, 0 ]
end
function F1(x)
    r, v, m = x
    return [ 0, Tmax/m, -b*Tmax ]
end

# Parameters for the dynamics
Cd = 310
β = 500
b = 2
Tmax_0 = 3.5
Tmax_f = 1.0

# Solve the problem with a reference value of Tmax
Tmax = Tmax_0
sol0 = solve(goddard; display=false)
@printf("Objective for reference solution: %.6f\n", objective(sol0))
```

Then, we perform the continuation on the maximal thrust.

```@example main-cont
sol = sol0 # Initialize the solution with the reference solution
data = DataFrame(Tmax=Float64[], Objective=Float64[], Iterations=Int[])
for Tmax_local ∈ range(Tmax_0, Tmax_f, length=6)
    global Tmax = Tmax_local # Update the global variable Tmax
    global sol = solve(goddard; init=sol, display=false)
    push!(data, (Tmax=Tmax, Objective=objective(sol), Iterations=iterations(sol)))
end 
println(data)
```

We plot now the objective with respect to the maximal thrust, as well as both solutions for `Tmax=3.5` and `Tmax=1`.

```@example main-cont
using Plots.PlotMeasures # for leftmargin

plt_obj = plot(data.Tmax, data.Objective;
    seriestype=:scatter,
    title="Goddard problem",
    label="r(tf)", 
    xlabel="Maximal thrust (Tmax)",
    ylabel="Maximal altitude r(tf)")

plt_sol = plot(sol0; label="Tmax="*string(data.Tmax[1]))
plot!(plt_sol, sol;  label="Tmax="*string(data.Tmax[end]))

layout = grid(2, 1, heights=[0.2, 0.8])
plot(plt_obj, plt_sol; layout=layout, size=(800, 1000), leftmargin=5mm)
```
