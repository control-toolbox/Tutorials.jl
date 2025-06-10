
# [Free Initial and Final Times](@id tutorial-free-times)

In this tutorial, we explore optimal control problems with free initial time `t₀` and/or final time `t_f`. These problems require special treatment in both direct and indirect methods, particularly when handling time-varying constraints and objectives. 


```@example free-time
using OptimalControl
using NLPModelsIpopt
using BenchmarkTools
using DataFrames
using Plots


function double_integrator_mintf()
    @def ocp begin
        tf ∈ R, variable
        t ∈ [0, tf], time
        x ∈ R², state
        u ∈ R, control
        -1 ≤ u(t) ≤ 1
        x(0) == [0, 0]
        x(tf) == [1, 0]
        0.05 ≤ tf ≤ Inf
        ẋ(t) == [x₂(t), u(t)]
        tf → min
    end
    return ocp
end
nothing # hide
```

##  Pour la résolution directe on obtient .


```@example main-disc
sol = solve(ocp; disc_method=:trapeze, display=false)
plot(sol; size=(800, 800))
```

