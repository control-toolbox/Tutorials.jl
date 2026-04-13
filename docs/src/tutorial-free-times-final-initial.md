# [Free initial and final times](@id tutorial-free-times-final-initial)

```@meta
Draft = false
```

This tutorial is part of a series on optimal control problems with free time variables.
See also: [Free final time](@ref tutorial-free-times-final) and [Free initial time](@ref tutorial-free-times-initial).

In this tutorial, we consider an optimal control problem with the initial and final times as free variables, that is there are parts of the variable to be optimized. The required packages for the tutorial are:

```@example free-both-times
using OptimalControl          # Main package
using NLPModelsIpopt          # Direct solver
using Plots                   # Visualization
```

The problem we consider is the following:

```@example free-both-times

@def ocp begin

    v = (t0, tf) ∈ R², variable # the initial and final times are free variables

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

We now solve the problem using a direct method.

```@example free-both-times
sol = solve(ocp; grid_size=100)
nothing # hide
```

```@example free-both-times
sol # hide
```

And plot the solution.

```@example free-both-times
plot(sol)
```
