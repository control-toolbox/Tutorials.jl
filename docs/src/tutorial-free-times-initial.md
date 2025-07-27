```@meta
Draft = false
```

# [Optimal control problem with free initial time](@id tutorial-free-times-initial)

In this tutorial, we explore an optimal control problem with free initial time `t0`.

Here are the required packages for the tutorial:

```@example initial_time
using LinearAlgebra: norm
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

    t0 → max
end
nothing # hide
```

#  Direct resolution

Let us solve the problem with a direct method.

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

## Indirect method 


Define the pseudo-Hamiltonian and switching function
```@example initial_time

H(x, p, u) = p[1]*x[2] + p[2]*u

# Hamiltonian lift for switching condition
F1(x) = [0, 1]
H1 = Lift(F1)

# Define flows for u = +1 and u = -1
const u_pos = 1
const u_neg = -1


f_pos = Flow(ocp, (x, p, t0) -> u_pos)
f_neg = Flow(ocp, (x, p, t0) -> u_neg)
```

Shooting function adapted for free initial time
```@example initial_time

function shoot!(s, p0, t1, t0)
    # Initial conditions
    x0 = [0, 0]  # fixed initial state
    tf = 0       # fixed final time
    
    # The flows (time runs from t0 to tf = 0)
    x1, p1 = f_pos(t0, x0, p0, t1)  # from t0 to t1
    xf, pf = f_neg(t1, x1, p1, tf)  # from t1 to tf = 0
    
    # Final control
    uf = -1
    
    # Target final state
    xf_target = [1, 0]
    
    # Shooting conditions
    s[1:2] = xf .- xf_target        # reach the target at tf = 0
    s[3]   = H(xf, pf, uf) - 1      # final condition on pseudo-Hamiltonian
    s[4]   = H1(x1, p1)             # switching condition at t1
end
```

Get initial guess from direct solution
```@example initial_time

t = time_grid(sol)
x = state(sol)
p = costate(sol)
φ(t) = H1(x(t), p(t))  # switching function

p0 = p(t[1])  # initial costate (at t0)
t1 = t[argmin(abs.(φ.(t)))]  # switching time
t0 = t[1]     # initial time (negative value)

println("p0 = ", p0)
println("t1 = ", t1)
println("t0 = ", t0)

#Test shooting function at initial guess
s = zeros(4)
shoot!(s, p0, t1, t0)
println("\nNorm of the shooting function: ‖s‖ = ", norm(s), "\n")

# Solve the nonlinear system
using NonlinearSolve

# Aggregated function
nle! = (s, ξ, λ) -> shoot!(s, ξ[1:2], ξ[3], ξ[4])

# Initial guess: [p0[1], p0[2], t1, t0]
ξ_guess = [p0..., t1, t0]
prob = NonlinearProblem(nle!, ξ_guess)

# Solve the problem
indirect_sol = solve(prob; show_trace=Val(true), abstol=1e-8, reltol=1e-8)
# Extract solution
p0_opt = indirect_sol.u[1:2]
t1_opt = indirect_sol.u[3]
t0_opt = indirect_sol.u[4]
```
We can now plot and compare with the direct method.# Compute and plot the optimal trajectory

```@example initial_time

x0 = [0, 0]
tf = 0

# Concatenate flows: u=+1 from t0 to t1, then u=-1 from t1 to tf
f = f_pos * (t1_opt, f_neg)
flow_sol = f((t0_opt, tf), x0, p0_opt; saveat=range(t0_opt, tf, 100))

# Plot comparison
plt = plot(sol; label="direct", size=(800, 800))
plot!(plt, flow_sol; label="indirect", color=:red)
```