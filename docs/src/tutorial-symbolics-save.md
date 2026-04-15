# Cart-Pole Limit Cycle via Symbolic Lagrangian Mechanics

```@meta
Draft = true
```

This tutorial demonstrates how to combine **symbolic derivation of equations of motion** (via
[Symbolics.jl](https://symbolics.juliasymbolics.org/)) with **direct optimal control**
(via [OptimalControl.jl](https://control-toolbox.org/OptimalControl.jl/)) to find a
**periodic orbit** (limit cycle) for a cart-pole system.

The key idea is to let the computer do the hard mechanics: we write the Lagrangian in a few
lines, and the Euler–Lagrange equations — including mass-matrix inversion — are derived
automatically.

## The Cart-Pole System

The system consists of a cart of mass ``m_c`` sliding on a frictionless horizontal rail, with
a rigid pendulum of mass ``m_p`` and length ``l`` attached to it. A horizontal force ``u``
(the control input) acts on the cart.

```@raw html
<figure>
  <img src="cartpole_sketch.svg" alt="Cart-pole diagram" width="350"/>
  <figcaption>Fig. 1 — Cart-pole system. The angle θ is measured from the upright position.</figcaption>
</figure>
```

The configuration vector is ``q = (x,\, \theta)^\top``, where ``x`` is the cart position and
``\theta = 0`` corresponds to the **upright** (unstable) equilibrium of the pendulum.

### Positions

The Cartesian positions of the two bodies are:

```math
p_c = \begin{pmatrix} x \\ 0 \end{pmatrix}, \qquad
p_p = \begin{pmatrix} x + l\sin\theta \\ l\cos\theta \end{pmatrix}.
```

### Lagrangian

The kinetic and potential energies are:

```math
T = \tfrac{1}{2}m_c\,\|\dot{p}_c\|^2 + \tfrac{1}{2}m_p\,\|\dot{p}_p\|^2
  = \tfrac{1}{2}(m_c+m_p)\dot{x}^2
    + m_p l\,\dot{x}\dot{\theta}\cos\theta
    + \tfrac{1}{2}m_p l^2\dot{\theta}^2,
```

```math
V = m_p\,g\,l\cos\theta.
```

The Lagrangian is ``\mathcal{L} = T - V``, and the virtual work of the control force gives the
generalised force vector ``Q = (u,\, 0)^\top``.

### Euler–Lagrange Equations

The equations of motion follow from:

```math
\frac{d}{dt}\frac{\partial \mathcal{L}}{\partial \dot{q}_i}
- \frac{\partial \mathcal{L}}{\partial q_i} = Q_i, \qquad i = 1,2.
```

They can be written in the standard **manipulator form**:

```math
M(q)\,\ddot{q} = -C(q,\dot{q}) + \tau(q, u),
```

where the symmetric positive-definite **mass matrix** is:

```math
M(q) =
\begin{pmatrix}
  m_c + m_p & m_p l \cos\theta \\
  m_p l \cos\theta & m_p l^2
\end{pmatrix},
```

and the right-hand side collects Coriolis/gravity terms and the control torque. Instead of
deriving these by hand, we rely on Symbolics.jl to compute and **analytically invert**
``M(q)`` for us.

### State-Space Form

Defining the state ``X = (x,\,\theta,\,\dot{x},\,\dot{\theta})^\top``, the equations of
motion become the first-order system:

```math
\dot{X}(t) = f\!\left(X(t),\, u(t)\right) =
\begin{pmatrix}
  \dot{x} \\ \dot{\theta} \\ M^{-1}(q)\bigl(-C(q,\dot{q}) + \tau(q,u)\bigr)
\end{pmatrix}.
```

## The Optimal Control Problem

We look for a **limit cycle** of the nonlinear dynamics: a trajectory that returns exactly to
its initial condition after a fixed period ``t_f``. The cost penalises the total control
energy:

```math
\min_{u(\cdot)}\; \int_0^{t_f} u(t)^2\,\mathrm{d}t
```

```math
\text{subject to} \quad \dot{X}(t) = f(X(t), u(t)), \quad t \in [0, t_f],
```

```math
X(t_f) = X(0).
```

The periodicity constraint ``X(t_f) = X(0)`` — combined with a non-trivial initial condition
— forces the solver to find an orbit rather than the trivial rest solution.

## Implementation

### Setup & Imports

```@example main
# ==========================================
# Setup & Imports
# ==========================================
import Pkg
Pkg.activate("Control.jl/ModelingToolkit.jl")

using OptimalControl
using Plots
using StaticArrays

import NLPModelsIpopt

using Symbolics
```

### Physical Parameters and Symbolic Variables

We declare all parameters both as numerical constants (for the final function evaluation) and
as symbolic variables (for the Lagrangian computation).

```@example main
# ==========================================
# 1. Variables and Parameters
# ==========================================
# Physical Constants
const m_c_val = 5.0
const m_p_val = 1.0
const l_val = 2.0
const g_val = 9.81
const tf_val = 2.0

# Symbolic Variables
@variables t
D = Differential(t)
@variables m_c m_p l g u
@variables x(t) θ(t)
@variables v ω dv dω # Static variables to solve for accelerations

q = [x, θ]
nothing # hide
```

### Automated Kinematics and Lagrangian

We express the positions, kinetic energy, potential energy, and virtual work symbolically.
The time derivatives ``\dot{p}_c``, ``\dot{p}_p`` are computed automatically by `D.(...)`.

```@example main
# ==========================================
# 2. Automated Kinematics & Jacobians
# ==========================================
p_c = [x, 0.0]
p_p = [x + l * sin(θ), l * cos(θ)]
F = [u, 0.0]

T = 0.5 * m_c * sum(D.(p_c) .^ 2) + 0.5 * m_p * sum(D.(p_p) .^ 2)
V = g * (m_p * p_p[2])
P_non_conservative = transpose(D.(p_c)) * F
nothing # hide
```

### Euler–Lagrange Equations and Mass-Matrix Inversion

Starting from ``\mathcal{L} = T - V``, Symbolics.jl computes the three terms of the
Euler–Lagrange equations and assembles the residual vector. Substituting static aliases
``(v,\omega,\dot{v},\dot\omega)`` for the time derivatives makes it possible to identify
the mass matrix ``M`` as the Jacobian of the residual with respect to the accelerations
``(\dot{v}, \dot\omega)``. The system ``M\,a = -b`` is then solved analytically.

```@example main
# ==========================================
# 3. Euler-Lagrange & Matrix Inversion
# ==========================================
L = T - V

A = D.(Symbolics.gradient(L, D.(q)))
B = Symbolics.gradient(L, q)
Q = Symbolics.gradient(P_non_conservative, D.(q))

# Euler-Lagrange: d/dt( dL/dq_dot ) - dL/dq = Forces
el_eqs = expand_derivatives.(A - B - Q)

# Freeze time derivatives into static algebraic variables
sub_rules = Dict(D(x) => v, D(θ) => ω, D(D(x)) => dv, D(D(θ)) => dω)
res = Symbolics.substitute.(el_eqs, (sub_rules,))

# Extract Mass Matrix and Bias vector
Mass = Symbolics.jacobian(res, [dv, dω])
bias = Symbolics.substitute.(res, (Dict(dv => 0.0, dω => 0.0),))

# Analytically invert the Mass Matrix (Symbolics handles 2x2 natively)
accel = Mass \ (-bias)
accel = Symbolics.simplify_fractions.(accel)

# The fully explicit state derivatives: X_dot = [v, ω, v_dot, ω_dot]
dx_dt = [v, ω, accel[1], accel[2]]
nothing # hide
```

### Code Generation

`build_function` compiles the symbolic expression `dx_dt` into a native Julia function.
The `force_SA=true` flag generates a **StaticArrays** kernel, which avoids heap allocations
inside the ODE right-hand side — important for solver performance.

```@example main
# ==========================================
# 4. Extract Explicit Dynamics into a Julia Function
# ==========================================
f_expr = build_function(dx_dt, [x, θ, v, ω], u, [m_c, m_p, l, g];
    expression=Val{false}, force_SA=true)
f_mtk = f_expr[1]

const p_vals = [m_c_val, m_p_val, l_val, g_val]

cartpole_dynamics(X, U) = f_mtk(X, U, p_vals)
nothing # hide
```

### Optimal Control Problem Definition

We now formulate the optimal control problem using the `@def` macro from OptimalControl.jl.
The periodicity condition is encoded as ``X(t_f) - X(0) = 0``, and the dynamics plug
directly into the `Ẋ(t) == ...` constraint.

```@example main
# ==========================================
# 5. OptimalControl.jl dynamic optimization problem
# ==========================================
@def cartpole_ocp begin
    t ∈ [0, tf_val], time
    X = (x_ocp, θ_ocp, v_ocp, ω_ocp) ∈ R⁴, state
    F_ctrl ∈ R, control

    X(0) == [0, 0, 0, 0.1]
    X(tf_val) - X(0) == [0, 0, 0, 0]

    # Inject our explicitly inverted Symbolic dynamics!
    Ẋ(t) == cartpole_dynamics(X(t), F_ctrl(t))

    ∫(F_ctrl(t)^2) → min
end
```

### Solving the NLP

The problem is transcribed into a nonlinear program using direct collocation on a uniform
grid of 100 intervals, then handed to **Ipopt** via NLPModelsIpopt.

```@example main
# ==========================================
# 6. Solvers and Collocation
# ==========================================
initial_guess = (state=[0.0, 0.0, 0.0, 0.1], control=0.0)

sol = solve(cartpole_ocp; display=true, grid_size=100, init=initial_guess)
```

### Results

```@example main
# ==========================================
# 7. Extracting Results
# ==========================================
println("--- Optimal Limit Cycle Found ---")
println("Total Control Energy (∫F² dt): ", sol.objective)

tsol = time_grid(sol)
# tsol = range(sol.model.times.initial.time, sol.model.times.final.time, 1001)
Xsol = state(sol).(tsol)
Fsol = control(sol).(tsol)

xsol = [X[1] for X in Xsol]
θsol = [X[2] for X in Xsol]
vsol = [X[3] for X in Xsol]
ωsol = [X[4] for X in Xsol]

qplot = plot(tsol, [xsol θsol], label=["x" "θ"], title="Configuration")
dqplot = plot(tsol, [vsol ωsol], label=["v" "ω"], title="Velocities")
uplot = plot(tsol, Fsol, label="u", title="Control", linetype=:steppost)

plot(qplot, dqplot, uplot, layout=3, size=(800, 600))
```

The three panels show, from top to bottom: the cart position ``x`` and pendulum angle
``\theta``, the corresponding velocities ``\dot{x}`` and ``\dot\theta``, and the optimal
control force ``u``. The periodicity of the state trajectory confirms that a genuine limit
cycle has been found.
