```@meta
Draft = false
```

# [Free Initial and Final Times](@id tutorial-free-times)

In this tutorial, we explore optimal control problems with free initial time `t₀` and/or final time `t_f`. These problems require special treatment in both direct and indirect methods, particularly when handling time-varying constraints and objectives.


```@example main-disc
using OptimalControl
using NLPModelsIpopt
using BenchmarkTools
using DataFrames
using Plots
using Printf

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

To allow the use of a free final time, you define t0 and/or tf as variables, and then set the global time t between these two variables. Here is a complete example :

```julia
v=(t0, tf) ∈ R², variable
t ∈ [t0, tf], time
```

You may need to add this type of constrain in your problem definition :

```julia
0.05 ≤ tf ≤ Inf
```

It is to ensure and help the convergence of the algorithm depending on your problem and it can also be needed when dealing with problems having a free initial time, as we will see in the next example.

##  Direct resolution with free final time :

We now solve the problem using a direct method, with automatic treatment of the free final time.

```@example main-disc
ocp = double_integrator_mintf()
sol = solve(ocp; grid_size=100)
plot(sol; label="direct", size=(800, 800))
```

## Verification of results

Here the theorical part :
```math
H = p1*x1 + p2*u
```
Conditions of Pontryagin's theorem :
```math
p1' = 0, \quad p2' = -p1 \\
p1 = p1(0) = constante = 1\\
p2 = -p1*t + p2(0) = 1-t\\
```
Transversal conditions :
```math
H[tf] = -p° = 1
```
We find u :
```math
u(t) = 1 ∈ [0,t1] \\
u(t) = -1 ∈ [t1, tf]
```
Where t1 is a constant to determined.

Now we have to integrate x' :

On t ∈ [0,t1] :
```math
x1' = x2 \\
x2'= u = 1 \\
```

```math
x2(t) = t \\
x1(t) = (1/2)*t^2
```

When t = t1 :
```math
x2(t1) = t1 \\
x1(t1) = (1/2)*t1^2
```

On t ∈ [t1,tf]
```math
x2(t) = -t + 2*t1 \\
x1(t) = -(1/2)*t^2 + 2*t1*t + C2 \\
```

```math
x1(t1) = -(1/2)*t1^2 + 2*t1^2 + C2 = (1/2)*t1^2 \\
C2 = -t1^2 \\
```

```math
x1(t) = -(1/2)*t^2 + 2*t1*t - t1^2
```

Finally you can solve the terminal conditions :
```math
x1(tf) = 1, \quad x2(tf) = 0 \\
```

```math
x2(tf) = -tf + 2*t1 = 0 \\
tf = 2*t1 \\
```

```math
x1(tf) = -(1/2)*tf^2 + 2*t1*tf - t1^2 = 1 \\
-2*t1^2 + 4*t1 - t1^2 = t1 = 1 \\
```

```math
t1 = 1, \quad t2 = 2
```

To sum up we find the following solutions :
```math
x1(t) = \begin{cases}
(1/2)*t^2 & \text{si } t \in [0, 1) \\
-(1/2)*t^2 + 2*t - 1 & \text{si } t \in [1, 2]
\end{cases}
\qquad
x2(t) = \begin{cases}
t & \text{si } t \in [0, 1) \\
2 - t & \text{si } t \in [1, 2]
\end{cases}
```

```math
p1(t) = 1, \quad p2(t) = 1-t, \quad p°=-1
```

```math
u(t) = \begin{cases}
1 & \text{si } t \in [0, 1) \\
-1 & \text{si } t \in [1, 2]
\end{cases}
```

We can analyse the influence of using different discretization sizes (grid_size), and observed the following results for the optimal tf:

```@example main-disc
for N in [20, 50, 100, 200]
    solN = solve(ocp; grid_size=N, display=false)
    @printf("grid_size = %3d → tf = %.5f\n", N, objective(solN))
end
```

This example shows that problems with a free final time can be sensitive to discretization. A small grid may lead to suboptimal or slightly inaccurate results.


# Now we will try an example with a free initial time

```@example initial_time
using OptimalControl
using NLPModelsIpopt
using BenchmarkTools
using DataFrames
using Plots
using Printf

function double_integrator_mint0()
    @def ocp begin
        t0 ∈ R, variable
        tf=0
        t ∈ [t0, tf], time
        x ∈ R², state
        u ∈ R, control
        -1 ≤ u(t) ≤ 1
        x(t0) == [0, 0]
        x(tf) == [1, 0]
        0.05 ≤ -t0 ≤ Inf
        ẋ(t) == [x₂(t), u(t)]
        -t0 → min
    end
    return ocp
end
nothing # hide
```

#  Direct resolution with free initial time :

We now solve the problem using a direct method, with automatic treatment of the free initial time.

```@example initial_time
ocp = double_integrator_mint0()
sol = solve(ocp; grid_size=100)
plot(sol; label="direct", size=(800, 800))
```








## From free final time to fixed final time :
