```@meta
Draft = false
```

# [Optimal control problem with free initial and free final times](@id tutorial-free-times-final-initial)

In this tutorial, we explore optimal control problems with free initial time `t₀` and final time `t_f`. 


## Here are the required packages for the tutorial:

```@example both_time
using OptimalControl
using NLPModelsIpopt
using BenchmarkTools
using DataFrames
using Plots
using Printf
```
## Definition of the problem

```@example both_time

@def ocp begin
    v=(t0, tf) ∈ R², variable
    t ∈ [t0, tf], time
    x ∈ R², state
    u ∈ R, control
    -1 ≤ u(t) ≤ 1
    x(t0) == [0, 0]
    x(tf) == [1, 0]
    0.05 ≤ t0 ≤ 10
    0.05 ≤ tf ≤ 10
    0.01 ≤ tf - t0 ≤ Inf
    ẋ(t) == [x₂(t), u(t)]
    t0 → max
end

nothing # hide
```

#  Direct resolution with both final and inital times being free:


We now solve the problem using a direct method, with automatic treatment of the free initial time.

```@example both_time
sol = solve(ocp; grid_size=100)
```
And plot the solution.

```@example both_time
plot(sol; label="direct", size=(800, 800))
```
