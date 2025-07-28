```@meta
Draft = false
```

# [Optimal control problem with free initial time](@id tutorial-free-times-initial)

In this tutorial, we explore an optimal control problem with free initial time `t0`.

## Definition of the problem

We consider the double integrator in minimum time. We start from the initial condition $x(t_0) = [0, 0]$ and we aim to reach at the final time $t_f = 0$ the target $x(t_f) = [1, 0]$. The final time is fixed but the initial time is free. The objective is thus to maximise `t0` which is negative.

```@example initial_time
using OptimalControl

tf = 0 # final time

@def ocp begin

    t0 ∈ R, variable        # the initial time is free
    t ∈ [t0, tf], time
    x ∈ R², state
    u ∈ R, control

    -1 ≤ u(t) ≤ 1

    x(t0) == [0, 0]
    x(tf) == [1, 0]

    0.05 ≤ -t0 ≤ Inf        # t0 is negative

    ẋ(t) == [x₂(t), u(t)]

    -t0 → min               # or t0 → max

end
nothing # hide
```

#  Direct resolution

Let us solve the problem with a direct method.

```@example initial_time
using NLPModelsIpopt
sol = solve(ocp)
nothing # hide
```

The initial time obtained with the direct method is:

```@example initial_time
t0 = variable(sol)
```

We can plot the solution.

```@example initial_time
using Plots
plot(sol; label="direct", size=(700, 700))
```

## Mathematical verification

Let us compute the solution analytically, thanks to the Pontryagin's Maximum Principle (PMP). First, we define the pseudo-Hamiltonian

```math
H(x, p, u) = p_1 x_2 + p_2 u.
```

If $(x(\cdot), u(\cdot), t_0)$ is a solution, then, there exists $(p(\cdot), p^0) \ne 0$, $p^0 \le 0$, satisfying the conditions given by the (PMP). First, since the initial time is free, we have $H(x(t_0), p(t_0), u(t_0)) = -p^0$. Then, the costate $p(\cdot)$ satisfies the adjoint equations

```math
    \begin{aligned}
        \dot{p}_1(t) &= 0, \\
        \dot{p}_2(t) &= -p_1(t),
    \end{aligned}
```

which gives $p_1(t) = \alpha$, $p_2(t) = -\alpha t + \beta$. The maximisation condition from the PMP implies that $u(t) = 1$ if $p_2(t)>0$ and $u(t)=-1$, if $p_2(t)<0$. We switch at time $t_s$ from one control to the other if $p_2(t_s) = 0$, that is when $\beta = \alpha\, t_s$. Note that $p_2$ cannot vanish on an interval of non-empty interior since otherwise, we would have $\alpha = \beta = 0$ and so $H(x(t_0), p(t_0), u(t_0)) = -p^0 = 0$ but $(p(\cdot), p^0) \ne 0$. Hence, there is at most one switching time $t_s$. We can demonstrate that in the setting above, there is exactly one switching time and we switch from $u=-1$ to $u=+1$. Thus, the optimal control is of the form $u(t) =  1$ for $t \in [t_0, t_s)$ and $u(t) = -1$ for $t \in [t_s, 0]$.

Let us find now $\alpha$, $\beta$, $t_s$ and $t_0$. We can integrate the system: On the interval $[t_0, t_s]$ we have $\dot{x}_2(t) = u(t) = 1$, hence, $x_2(t) = t - t_0$. Since $\dot{x}_1(t) = x_2(t)$, we get 

```math
x_1(t) = \frac{(t - t_0)^2}{2}.
```

At switching time $t_s$ whe have:

```math
x_2(t_s) = t_s - t_0, \quad x_1(t_s) = \frac{(t_s - t_0)^2}{2},
```

that are the initial conditions of the integrations on the interval $[t_s, t_f]$:

```math
\dot{x}_2(t) = u(t) = -1 \quad \Rightarrow \quad x_2(t) = x_2(t_s) - (t - t_s) = 2\, t_s - t_0 - t
```

and

```math
\dot{x}_1(t) = x_2(t) \quad \Rightarrow \quad x_1(t) = x_1(t_s) + \int_{t_s}^t x_2(\sigma)\, \mathrm{d}\sigma = 
\frac{(t_s - t_0)^2}{2} + (2\, t_s - t_0) (t - t_s) - \frac{t^2}{2} + \frac{t_s^2}{2}.
```

The final condition on $x_2$ gives ($t_f = 0$):

```math
x_2(0) = 0 \quad \Rightarrow \quad t_0 = 2\, t_s.
```

The final condition of $x_1$ gives:

```math
x_1(0) = t_s^2 = 1 \quad \Rightarrow \quad t_s = -1.
```

We deduce:

```math
t_0 = 2\, t_s = -2.
```

To get $\alpha$ and $\beta$ we use the initial condition on the pseudo-Hamiltonian and the swtiching condition $p_2(t_s) = 0$. The switching condition gives $\beta = \alpha\, t_s = -\alpha$. The initial condition on the pseudo-Hamiltonian gives:

```math
H(x(t_0), p(t_0), u(t_0)) = p_1(t_0) x_2(t_0) + p_2(t_0) u(t_0) = -\beta = -p^0 = 1.
```

The dual variable $p^0 \ne 0$ since otherwise $(p(\cdot), p^0) = 0$ which is not. Then, $p^0 < 0$ and we can normalise it to $p^0 = -1$. At the end, we have:

```math
\alpha = 1, \quad \beta = -1, \quad t_s = -1, \quad t_0 = -2.
```

## Indirect method 

Let us compute numerically the solution of the PMP using the solution from the direct method. We define the pseudo-Hamiltonian and the flows associated to the two controls.

```@example initial_time
# pseudo-Hamiltonian
H(x, p, u) = p[1]*x[2] + p[2]*u

# Define flows for u = +1 and u = -1
const u_pos = 1
const u_neg = -1

using OrdinaryDiffEq
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
using LinearAlgebra: norm
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