# [Optimal control problem with free initial time](@id tutorial-free-times-initial)

In this tutorial, we study a **minimum-time optimal control problem with free initial time** $t_0$ (negative). The goal is to determine the latest possible starting time so that the system reaches a fixed final state at $t_f = 0$.

The system is the classic **double integrator**:

- Initial state: $x(t_0)=[0,0]$, final state: $x(t_f)=[1,0]$  
- Control bounds: $u(t) \in [-1,1]$  
- Objective: maximize $t_0$ (equivalently minimize $-t_0$)  
- Dynamics: $\dot{x}_1 = x_2, \ \dot{x}_2 = u$

```@example initial_time
using LinearAlgebra: norm
using OptimalControl
using NLPModelsIpopt
using NonlinearSolve
using OrdinaryDiffEq
using Plots
using Printf
```

## Problem definition

We consider the following setup:

```@example initial_time
tf = 0          # Fixed final time
x0 = [0, 0]     # Initial state
xf = [1, 0]     # Final state

@def ocp begin

    t0 ∈ R, variable           # Free initial time
    t ∈ [t0, tf], time
    x ∈ R², state
    u ∈ R, control

    # Control bounds
    -1 ≤ u(t) ≤ 1

    # Boundary conditions
    x(t0) == x0
    x(tf) == xf

    # Ensure t0 is negative and bounded for numerical stability
    0.05 ≤ -t0 ≤ Inf

    # Dynamics
    ẋ(t) == [x₂(t), u(t)]

    # Objective: maximize t0 (equivalent to minimizing -t0)
    -t0 → min

end
nothing # hide
```

This explicitly declares the **free initial time** via `t0 ∈ R, variable` and defines the **time interval** `[t0, tf]`.

## Direct solution

To improve convergence of the direct solver, we constrain `t0` as follows:

```julia
-Inf ≤ t0 ≤ -0.05
```

We solve the problem using a **direct transcription method**:

```@example initial_time
sol = solve(ocp; grid_size=100)
nothing # hide
```

The solution can be visualized:

```@example initial_time
plt = plot(sol; label="direct", size=(800, 600))
```

## Verification of the results

Using **Pontryagin's Maximum Principle (PMP)**, the optimal control is **bang-bang** with a single switch at $t_1$:

```math
u(t) =
\begin{cases}
1, & t \in [t_0,t_1),\\
-1, & t \in [t_1,0].
\end{cases}
```

The corresponding trajectory, costate, switching time, and initial time are:

```math
x_1(t) =
\begin{cases}
\tfrac{1}{2}(t-t_0)^2, & t \in [t_0,-1),\\
-\tfrac{1}{2}t^2 + 2 t_1 t - t_1^2 - t_0 t_1, & t \in [-1,0],
\end{cases}
\quad
x_2(t) =
\begin{cases}
t - t_0, & t \in [t_0,-1),\\
2 t_1 - t_0 - t, & t \in [-1,0],
\end{cases}
```

```math
p(t_0) = (1,1), \quad p(t_f) = (1,-1), \quad t_1=-1, \quad t_0=-2, \quad p^0=-1.
```

We can now compare the direct numerical solution with this theoretical result:

```@example initial_time
t0 = variable(sol)
u = control(sol)
p = costate(sol)
x = state(sol)
H(t) = p(t)[1]*x(t)[2] + p(t)[2]*u(t)

@printf("t0 = %.5f", t0); println(lpad("(expected -2)", 33))
@printf("H(t0) = %.5f", H(t0)); println(lpad("(expected 1)", 30))
@printf("H(tf) = %.5f", H(tf)); println(lpad("(expected 1)", 30))
@printf("x(t0) = [%.5f, %.5f]", x(t0)[1], x(t0)[2]); println(lpad("(expected [0, 0])", 23))
@printf("x(tf) = [%.5f, %.5f]", x(tf)[1], x(tf)[2]); println(lpad("(expected [1, 0])", 24))
@printf("p(t0) = [%.5f, %.5f]", p(t0)[1], p(t0)[2]); println(lpad("(expected [1, 1])", 24))
@printf("p(tf) = [%.5f, %.5f]", p(tf)[1], p(tf)[2]); println(lpad("(expected [1, -1])", 24))
```

The numerical results match the theoretical solution almost exactly.

## Indirect (shooting) solution

We now solve the PMP system numerically using an **indirect method** (shooting approach), using the direct solution as an initial guess.

```@example initial_time
# Pseudo-Hamiltonian
H(x, p, u) = p[1]*x[2] + p[2]*u

# Hamiltonian lift used to compute the switching function
F1(x) = [0, 1]
H1 = Lift(F1)
nothing # hide
```

Define the flows corresponding to the two control laws $u=+1$ and $u=-1$:

```@example initial_time
const u_pos = 1
const u_neg = -1

# Each flow depends on tf since it is part of the optimization variables
f_pos = Flow(ocp, (x, p, tf) -> u_pos)
f_neg = Flow(ocp, (x, p, tf) -> u_neg)
nothing # hide
```

The **shooting function** enforces:

- Continuity of state and costate  
- Satisfaction of boundary conditions  
- Switching condition ($p_2(t_1) = 0$)  
- Hamiltonian normalization at initial time

```@example initial_time
function shoot!(s, p0, t1, t0)
    x_t0 = x0
    p_t0 = p0

    # Forward under u = +1, then u = -1
    x_t1, p_t1 = f_pos(t0, x_t0, p_t0, t1)
    x_tf, p_tf = f_neg(t1, x_t1, p_t1, tf)

    u_t0 = 1  # Initial control

    # Shooting conditions:
    s[1:2] = x_tf - xf                  # Terminal state match
    s[3]   = H(x_t0, p_t0, u_t0) - 1    # Hamiltonian normalization
    s[4]   = H1(x_t1, p_t1)             # Switching condition φ(t₁)=0
end
nothing # hide
```

To help the nonlinear solver converge, we build a good initial guess from the direct solution:

```@example initial_time
t = time_grid(sol)
x = state(sol)
p = costate(sol)
φ(t) = H1(x(t), p(t))  # Switching function

p0 = p(t[1])
t1 = t[argmin(abs.(φ.(t)))]
t0 = t[1]

println("t0 = ", t0)
println("p0 = ", p0)
println("t1 = ", t1)

# Evaluate the norm of the shooting function at initial guess
s = zeros(4)
shoot!(s, p0, t1, t0)
println("\n‖s‖ (initial guess) = ", norm(s), "\n")
```

We can now solve the system using an **indirect shooting method**:

```@example initial_time
# Aggregated nonlinear system
shoot!(s, ξ, λ) = shoot!(s, ξ[1:2], ξ[3], ξ[4])

# Define the problem and initial guess
ξ_guess = [p0..., t1, t0]
prob = NonlinearProblem(shoot!, ξ_guess)

# Solve the nonlinear system
indirect_sol = solve(prob; show_trace=Val(true), abstol=1e-8, reltol=1e-8)
nothing # hide
```

```@example initial_time
indirect_sol # hide
```

Compare **indirect and direct solutions**:

```@example initial_time
p0 = indirect_sol.u[1:2]
t1 = indirect_sol.u[3]
t0 = indirect_sol.u[4]

# Reconstruct the full trajectory
f = f_pos * (t1, f_neg)
flow_sol = f((t0, tf), x0, p0; saveat=range(t0, tf, 200))

# Plot comparison
plot!(plt, flow_sol; label="indirect", color=:red, linestyle=:dash)
```
