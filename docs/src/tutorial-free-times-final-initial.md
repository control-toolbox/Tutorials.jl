# [Free initial and final times](@id tutorial-free-times-final-initial)

In this tutorial, we consider an optimal control problem with the initial and final times as free variables, that is there are parts of the variable to be optimized. The required packages for the tutorial are:

```@example both_time
using OptimalControl
using NLPModelsIpopt
using Plots
```

The problem we consider is the following:

```@example both_time

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

```@example both_time
sol = solve(ocp; grid_size=100)
nothing # hide
```

```@example both_time
sol # hide
```

And plot the solution.

```@example both_time
plot(sol)
```
