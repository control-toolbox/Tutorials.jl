```@meta
Draft = true
```

# [Optimal control problem with free initial time](@id tutorial-free-times-initial)

In this tutorial, we explore an optimal control problem with free initial time `t0`.

Here are the required packages for the tutorial:

```@example initial_time
using OptimalControl
using NLPModelsIpopt
using BenchmarkTools
using DataFrames
using Plots
using Printf
```

## Definition of the problem

We consider the double integrator in minimum time but we fix the final time and let the initial time free. The objective is thus to maximise `t0`.

```@example initial_time
tf = 0 # final time
@def ocp begin

    t0 ∈ R, variable    # the initial time is free
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
nothing # hide
```

#  Direct resolution

We now solve the problem using a direct method, with automatic treatment of the free initial time.

```@example initial_time
sol = solve(ocp)
```

The initial time solution is:

```@example initial_time
t0 = variable(sol)
```

We can plot the solution.

```@example initial_time
plot(sol; label="direct", size=(800, 800))
```

## Mathematical verification

Here is the theoretical part using Pontryagin's Maximum Principle:
```math
H = p_1 x_2 + p_2 u + 1
```

Conditions from Pontryagin’s theorem:
```math
    \begin{aligned}

        \dot{p}_1(t) &= 0,  \quad \Rightarrow \quad p_1 = c_1 \quad (\text{constant}) \\

        \dot{p}_2(t) &= -p_1 \quad \Rightarrow \quad p_2 = -c_1 t + c_2
    \end{aligned}
```

Switching condition:
```math
p_2(t_s) = 0 \quad \Rightarrow \quad c_2 = c_1 t_s
```

Optimal control:
```math
u(t) = 1 \quad \text{on} \quad [t_0, t_s] \\
u(t) = -1 \quad \text{on} \quad [t_s, 0]
```

Now we integrate the system:

On ( t in [t_0, t_s] ) :
```math
x_2' = u = 1 \quad \Rightarrow \quad x_2(t) = t - t_0 \\
x_1' = x_2 \quad \Rightarrow \quad x_1(t) = \frac{(t - t_0)^2}{2}
```

At switching time ( t = t_s ) :
```math
\dot{x}_2(t_s) = t_s - t_0 \\
\dot{x}_1(t_s) = \frac{(t_s - t_0)^2}{2}
```

when ( t in [t_s, 0] ) :
```math
\dot{x}_2(t) = u(t) = -1 \quad \Rightarrow \quad x_2(t) = x_2(t_s) - (t - t_s) \\
\dot{x}_1(t) = x_2(t) \quad \Rightarrow \quad x_1(t) = x_1(t_s) + \int_{t_s}^t x_2(s) ds
```

Final velocity condition:
```math
x_2(0) = 0 \quad \Rightarrow \quad t_s - t_0 + t_s = 0 \quad \Rightarrow \quad t_0 = 2 t_s
```

Final position:
```math
x_1(0) = x_1(t_s) + \frac{t_s^2}{2} \quad \Rightarrow \quad x_1(0) = t_s^2 = 1 \quad \Rightarrow \quad t_s = -1
```

We deduce:
```math
t_0 = 2 * t_s = -2
```

### Final solution:
- Switching time:
```math
t_s = -1 
```
- Initial time:
```math
t_0 = -2 
```

Control:
```math
u(t) = 1 \quad \text{on} \quad [-2, -1] \\
u(t) = -1 \quad \text{on} \quad [-1, 0]
```
