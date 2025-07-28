```@meta
Draft = false
```

# [Optimal control problem with free initial time](@id tutorial-free-times-initial)

In this tutorial, we explore a minimum-time optimal control problem where the **initial time $t_0$ is free** (but negative). We aim to determine the latest possible starting time for the system to still reach the goal at a fixed final time $t_f = 0$.

## Problem definition

We consider the classic **double integrator** system. The goal is to move from rest at position $0$ to position $1$ (also at rest), using bounded control:

- Initial condition: $x(t_0) = [0, 0]$,
- Final condition: $x(t_f) = [1, 0]$ with $t_f = 0$ fixed,
- Control bounds: $u(t) \in [-1, 1]$,
- The objective is to **maximize $t_0$**, i.e., minimize $-t_0$, which corresponds to a **minimum-time formulation** with fixed arrival time and free departure time.

```@example initial_time
using OptimalControl

tf = 0          # Final time
x0 = [0, 0]     # Initial state
xf = [1, 0]     # Final state

@def ocp begin

    t0 ∈ R, variable                # t₀ is free and will be optimized
    t ∈ [t0, tf], time
    x ∈ R², state
    u ∈ R, control

    -1 ≤ u(t) ≤ 1                   # Control bounds

    x(t0) == x0                     # Initial constraint
    x(tf) == xf                     # Final constraint

    0.05 ≤ -t0 ≤ Inf                # Ensure t₀ is negative and bounded

    ẋ(t) == [x₂(t), u(t)]           # Dynamics

    -t0 → min                       # Equivalent to maximizing t₀

end
nothing # hide
```

## Direct solution

We solve the problem using a **direct transcription method**.

```@example initial_time
using NLPModelsIpopt
direct_sol = solve(ocp)
nothing # hide
```

The optimal initial time obtained is:

```@example initial_time
t0 = variable(direct_sol)
```

We can now plot the solution:

```@example initial_time
using Plots
plt = plot(direct_sol; label="direct", size=(700, 700))
```

## Pontryagin's Maximum Principle

The optimal control appears to be **bang-bang**, switching once from $u = 1$ to $u = -1$. This motivates us to apply the Pontryagin's Maximum Principle (PMP) analytically to verify and interpret the solution.

We define the **pseudo-Hamiltonian**:

```math
H(x, p, u) = p_1 x_2 + p_2 u.
```

Let $(x(\cdot), u(\cdot), t_0)$ be the solution. By PMP, there exists $(p(\cdot), p^0) \neq 0$ with $p^0 \leq 0$ and $p(\cdot)$ absolutely continuous, such that $H(x(t_0), p(t_0), u(t_0)) = -p^0$. The adjoint equations are:

```math
\dot{p}_1(t) = 0, \quad \dot{p}_2(t) = -p_1(t),
```

so $p_1(t) = \alpha$, $p_2(t) = -\alpha (t-t_0) + \beta$. The maximization condition gives: 

```math
    \left\{
    \begin{array}{lll}
        u(t) =  1 & \text{ if } & p_2(t) > 0, \\[0.5em]
        u(t) = -1 & \text{ if } & p_2(t) < 0, \\[0.5em]
        u(t) \in  [-1, 1] & \text{ if } & p_2(t) = 0.
    \end{array}
    \right.
```

We assume one switch at $t_1$ with $p_2(t_1) = 0 \Leftrightarrow \beta = \alpha (t_1-t_0)$.

!!! note
    If $p_2(t) = 0$ on a nontrivial interval, then $p \equiv 0$ and $p^0 = 0$, violating PMP.

We integrate the system in two phases:
- On $[t_0, t_1]$ with $u = 1$:  
  $x_2(t) = t - t_0$ and $x_1(t) = \frac{1}{2}(t - t_0)^2$.
- At switch $t = t_1$:
  ```math
  x_2(t_1) = t_1 - t_0, \quad x_1(t_1) = \frac{(t_1 - t_0)^2}{2}
  ```
- On $[t_1, 0]$ with $u = -1$:  
  $x_2(t) = x_2(t_1) - (t - t_1) = 2 t_1 - t_0 - t$.
  Integrating gives:
  ```math
  x_1(t) = x_1(t_1) + (2t_1 - t_0)(t - t_1) - \frac{t^2}{2} + \frac{t_1^2}{2}
  ```

Final constraints:
```math
x_2(0) = 0 \Rightarrow t_0 = 2 t_1, \quad x_1(0) = 1 \Rightarrow t_1 = -1
```

Therefore: $t_0 = -2$.

To find $\alpha$ and $\beta$, use:
- Switching: $\beta = \alpha (t_1-t_0) = \alpha$,
- PMP: $H(x(t_0), p(t_0), u(t_0)) = -p^0 = \beta = 1$.

So the full solution is:
```math
p(t_0) = (1, 1), \quad t_1 = -1, \quad t_0 = -2.
```

## Indirect method

We now solve the PMP system **numerically** using the direct solution as an initial guess.

```@example initial_time
H(x, p, u) = p[1]*x[2] + p[2]*u

# Define flows for constant controls
u_pos = 1
u_neg = -1

using OrdinaryDiffEq
f_pos = Flow(ocp, (x, p, t0) -> u_pos)
f_neg = Flow(ocp, (x, p, t0) -> u_neg)
nothing # hide
```

The **shooting function** encodes:
- continuity of state and costate,
- satisfaction of boundary conditions,
- switching condition ($p_2(t_1) = 0$),
- pseudo-Hamiltonian condition at $t_0$.

```@example initial_time
function shoot!(s, p0, t1, t0)
    x1, p1 = f_pos(t0, x0, p0, t1, t0)
    x2, p2 = f_neg(t1, x1, p1, tf, t0)

    s[1:2] = x2 - xf                    # Final state match
    s[3]   = H(x0, p0, u_pos) - 1       # Hamiltonian normalization
    s[4]   = p1[2]                      # Switching condition
end
nothing # hide
```

Get an initial guess from the direct solution:

```@example initial_time
t = time_grid(direct_sol)
x = state(direct_sol)
p = costate(direct_sol)
φ(t) = p(t)[2]                     # Switching function

t0 = variable(direct_sol)
p0 = p(t0)                         # Initial costate
t1 = t[argmin(abs.(φ.(t)))]        # Approximate switching time

println("Initial guess:\n\np0 = ", p0, "\nt1 = ", t1, "\nt0 = ", t0)

s = zeros(4)
shoot!(s, p0, t1, t0)
using LinearAlgebra: norm
println("\n‖s‖ = ", norm(s))
```

Solve the shooting equations:

```@example initial_time
using NonlinearSolve

nle! = (s, ξ, λ) -> shoot!(s, ξ[1:2], ξ[3], ξ[4])
ξ_guess = [p0..., t1, t0]
prob = NonlinearProblem(nle!, ξ_guess)

nle_sol = solve(prob; show_trace=Val(true), abstol=1e-8, reltol=1e-8)

# Extract solution
p0 = nle_sol.u[1:2]
t1 = nle_sol.u[3]
t0 = nle_sol.u[4]

println("\nIndirect solution:\n\np0 = ", p0, "\nt1 = ", t1, "\nt0 = ", t0)

shoot!(s, p0, t1, t0)
println("\nFinal shooting residual: ‖s‖ = ", norm(s))
```

Compare with the direct solution:

```@example initial_time
f = f_pos * (t1, f_neg)     # Concatenate bang-bang flow
indirect_sol = f((t0, tf), x0, p0, t0; saveat=range(t0, tf, 200))
plot!(plt, indirect_sol; label="indirect", color=:red, linestyle=:dash)
```
