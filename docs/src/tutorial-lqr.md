# A simple Linear–quadratic regulator example

## Problem statement

We consider the following Linear Quadratic Regulator (LQR) problem, which consists in minimizing

```math
    \frac{1}{2} \int_{0}^{t_f} \left( x_1^2(t) + x_2^2(t) + u^2(t) \right) \, \mathrm{d}t 
```

subject to the dynamics

```math
    \dot x_1(t) = x_2(t), \quad \dot x_2(t) = -x_1(t) + u(t), \quad u(t) \in \R
```

and the initial condition

```math
    x(0) = (0,1).
```

We aim to solve this optimal control problem for different values of $t_f$.  

## Required packages

We begin by importing the necessary packages:

```@example main-lqr
using OptimalControl
using NLPModelsIpopt
using Plots
using Plots.PlotMeasures # for leftmargin, bottommargin
```

## Problem definition

We define the LQR problem as a function of the final time `tf`:

```@example main-lqr
function lqr(tf)

    x0 = [ 0
           1 ]

    ocp = @def begin
        t ∈ [0, tf], time
        x ∈ R², state
        u ∈ R, control
        x(0) == x0
        ẋ(t) == [x₂(t), - x₁(t) + u(t)]
        0.5∫( x₁(t)^2 + x₂(t)^2 + u(t)^2 ) → min
    end

    return ocp
end
nothing # hide
```

## Solving the problem for different final times

We solve the problem for $t_f \in \{3, 5, 30\}$.

```@example main-lqr
solutions = []   # empty list of solutions
tfs = [3, 5, 30]

for tf ∈ tfs
    solution = solve(lqr(tf), display=false)
    push!(solutions, solution)
end
nothing # hide
```

## Plotting the Solutions

We plot the state and control variables using normalized time $s = (t - t_0)/(t_f - t_0)$:

```@example main-lqr
plt = plot(solutions[1], :state, :control; time=:normalize, label="tf = $(tfs[1])")
for (tf, sol) ∈ zip(tfs[2:end], solutions[2:end])
    plot!(plt, sol, :state, :control; time=:normalize, label="tf = $tf")
end

px1 = plot(plt[1], legend=false, xlabel="s", ylabel="x₁")
px2 = plot(plt[2], legend=true, xlabel="s", ylabel="x₂")
pu  = plot(plt[3], legend=false, xlabel="s", ylabel="u")
plot(px1, px2, pu, layout=(1, 3), size=(800, 300), leftmargin=5mm, bottommargin=5mm)
```

!!! note "Nota bene"

    We can observe that $x(t_f)$ converges to the origin as $t_f$ increases.

## Known issue

The following definition will lead to an error when solving the problem. This is a [known issue](@extref OptimalControl manual-abstract-known-issues).

```@repl main-lqr

x0 = [ 0
       1 ]

A = [  0 1
      -1 0 ]

B = [ 0
      1 ]

Q = [ 1 0
      0 1 ]

R = 1

tf = 3

ocp = @def begin
    t ∈ [0, tf], time
    x ∈ R², state
    u ∈ R, control
    x(0) == x0
    ẋ(t) == A * x(t) + B * u(t)
    0.5∫( x(t)' * Q * x(t) + u(t)' * R * u(t) ) → min
end

solve(ocp)
```
