# [Optimal control problem with free final time](@id tutorial-free-times-final)

In this tutorial, we study an **optimal control problem with a free final time** `tf`. The following packages are required:

```@example main-free-final
using LinearAlgebra: norm
using NLPModelsIpopt
using NonlinearSolve
using OptimalControl
using OrdinaryDiffEq
using Plots
using Printf
```

## Definition of the problem

We consider the following setup:

```@example main-free-final
t0 = 0       # Initial time
x0 = [0, 0]  # Initial state
xf = [1, 0]  # Desired final state

@def ocp begin

    tf ∈ R, variable           # Final time is free (a decision variable)
    t ∈ [t0, tf], time         # Time interval depends on tf
    x ∈ R², state              # State vector x = [x₁, x₂]
    u ∈ R, control             # Scalar control

    # Control and time constraints
    -1 ≤ u(t) ≤ 1
    0.05 ≤ tf ≤ Inf            # Lower bound helps numerical convergence

    # Boundary conditions
    x(t0) == x0
    x(tf) == xf

    # Dynamics
    ẋ(t) == [x₂(t), u(t)]

    # Objective: minimize final time
    tf → min
end
nothing # hide
```

The **free final time** nature of the problem is explicitly stated by:

```julia
tf ∈ R, variable
```

and the time definition

```julia
t ∈ [t0, tf], time
```

indicates that the time horizon itself depends on the optimization variable `tf`.

## Direct resolution

To improve convergence of the direct solver, we constrain `tf` as follows:

```julia
0.05 ≤ tf ≤ Inf
```

We then solve the optimal control problem with a direct method, using an automatic handling of the free final time:

```@example main-free-final
sol = solve(ocp; grid_size=100)
nothing # hide
```

The solution can be visualized as:

```@example main-free-final
plt = plot(sol; label="direct", size=(800, 600))
```

## Verification of the results

Using **Pontryagin's Maximum Principle**, we know that the optimal trajectory consists of **two arcs**: one with $u=1$ and one with $u=-1$. Without further calculations, the optimal trajectory is:

```math
x_1(t) =
\begin{cases}
    \tfrac{1}{2}t^2, & t \in [0,1),\\
    -\tfrac{1}{2}t^2 + 2t - 1, & t \in [1,2],
\end{cases}
\quad
x_2(t) =
\begin{cases}
    t, & t \in [0,1),\\
    2-t, & t \in [1,2].
\end{cases}
```

The optimal control is:

```math
u(t) =
\begin{cases}
    1, & t \in [0,1),\\
-1, & t \in [1,2].
\end{cases}
```

And the costate, the switching time and the final time are:

```math
p_1(t)=1, \quad p_2(t)=1-t, \quad p^0=-1, \quad t_1 = 1, \quad t_f = 2.
```

We can now compare the direct numerical solution with this theoretical result:

```@example main-free-final
tf = variable(sol)
u = control(sol)
p = costate(sol)
x = state(sol)
H(t) = p(t)[1]*x(t)[2] + p(t)[2]*u(t)

@printf("tf = %.5f", tf); println(lpad("(expected 2)", 33))
@printf("H(t0) = %.5f", H(t0)); println(lpad("(expected 1)", 30))
@printf("H(tf) = %.5f", H(tf)); println(lpad("(expected 1)", 30))
@printf("x(t0) = [%.5f, %.5f]", x(t0)[1], x(t0)[2]); println(lpad("(expected [0, 0])", 23))
@printf("x(tf) = [%.5f, %.5f]", x(tf)[1], x(tf)[2]); println(lpad("(expected [1, 0])", 23))
@printf("p(t0) = [%.5f, %.5f]", p(t0)[1], p(t0)[2]); println(lpad("(expected [1, 1])", 24))
@printf("p(tf) = [%.5f, %.5f]", p(tf)[1], p(tf)[2]); println(lpad("(expected [1, -1])", 24))
```

The numerical results match the theoretical solution almost exactly.

## Indirect method

We now solve the same problem using an **indirect method** (shooting approach). We begin by defining the pseudo-Hamiltonian and the switching function.

```@example main-free-final
# Pseudo-Hamiltonian
H(x, p, u) = p[1]*x[2] + p[2]*u

# Hamiltonian lift used to compute the switching function
F1(x) = [0, 1]
H1 = Lift(F1)
nothing # hide
```

Define the flows corresponding to the two control laws $u=+1$ and $u=-1$:

```@example main-free-final
const u_pos = 1
const u_neg = -1

# Each flow depends on tf since it is part of the optimization variables
f_pos = Flow(ocp, (x, p, tf) -> u_pos)
f_neg = Flow(ocp, (x, p, tf) -> u_neg)
nothing # hide
```

The **shooting function** enforces the final and switching conditions:

```@example main-free-final
function shoot!(s, p0, t1, tf)
    x_t0 = x0
    p_t0 = p0

    # Forward integration under u = +1, then u = -1
    x_t1, p_t1 = f_pos(t0, x_t0, p_t0, t1)
    x_tf, p_tf = f_neg(t1, x_t1, p_t1, tf)

    u_tf = -1  # Final control

    # Shooting conditions:
    s[1:2] = x_tf - xf                # Terminal state must match xf
    s[3]   = H(x_tf, p_tf, u_tf) - 1  # Transversality condition
    s[4]   = H1(x_t1, p_t1)           # Switching condition φ(t₁)=0
end
nothing # hide
```

To help the nonlinear solver converge, we build a good initial guess from the direct solution:

```@example main-free-final
t = time_grid(sol)
x = state(sol)
p = costate(sol)
φ(t) = H1(x(t), p(t))  # Switching function

p0 = p(t0)
t1 = t[argmin(abs.(φ.(t)))]
tf = t[end]

println("p0 = ", p0)
println("t1 = ", t1)
println("tf = ", tf)

# Evaluate the norm of the shooting function at initial guess
s = zeros(4)
shoot!(s, p0, t1, tf)
println("\n‖s‖ (initial guess) = ", norm(s), "\n")
```

We can now solve the system using an **indirect shooting method**:

```@example main-free-final
# Aggregated nonlinear system
shoot!(s, ξ, λ) = shoot!(s, ξ[1:2], ξ[3], ξ[4])

# Define the problem and initial guess
ξ_guess = [p0..., t1, tf]
prob = NonlinearProblem(shoot!, ξ_guess)

# Solve the nonlinear system
indirect_sol = solve(prob; show_trace=Val(true), abstol=1e-8, reltol=1e-8)
nothing # hide
```

```@example main-free-final
indirect_sol # hide
```

Finally, we compare the indirect solution to the direct one:

```@example main-free-final
p0 = indirect_sol.u[1:2]
t1 = indirect_sol.u[3]
tf = indirect_sol.u[4]

# Reconstruct the full trajectory
f = f_pos * (t1, f_neg)
flow_sol = f((t0, tf), x0, p0; saveat=range(t0, tf, 200))

# Plot comparison
plot!(plt, flow_sol; label="indirect", color=2)
```
