# Cart-Pole Limit Cycle via Symbolic Lagrangian Mechanics

```@meta
Draft = false
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

The Lagrangian is ``\mathcal{L} = T - V``. The virtual work of the control force ``u``
acting on the cart gives the generalised force vector:

```math
W_\text{ctrl} = F \cdot \dot{p}_c = u\,\dot{x},
\quad\Longrightarrow\quad
Q = \frac{\partial W_\text{ctrl}}{\partial \dot{q}} = \begin{pmatrix} u \\ 0 \end{pmatrix}.
```

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
X(0) = X(t_f) = (0,\, 0,\, 0,\, 0.1)^\top.
```

The common boundary condition ``X(0) = X(t_f)`` — combined with a non-trivial initial
angular velocity — forces the solver to find an orbit rather than the trivial rest solution.

## Implementation

### Setup & Imports

```@setup main
# Project environment for this tutorial — edit the path for your local setup.
import Pkg
Pkg.activate(".")
```

```@example main
using OptimalControl
using Plots
using StaticArrays
using NLPModelsIpopt
using Symbolics
```

### Physical Parameters and Symbolic Variables

We declare all parameters both as numerical constants (for the final function evaluation) and
as symbolic variables (for the Lagrangian computation). Note that the symbolic time variable
`t` is used only inside the Symbolics.jl derivation and does not interfere with the time
variable `t` introduced later by the `@def` macro.

```@example main
# Physical constants
const m_c_val = 5.0
const m_p_val = 1.0
const l_val   = 2.0
const g_val   = 9.81
const tf_val  = 2.0

# Symbolic variables
@variables t
D = Differential(t)
@variables m_c m_p l g u
@variables x(t) θ(t)
@variables v ω dv dω   # static aliases for velocities and accelerations

q = [x, θ]
```

### Automated Kinematics and Lagrangian

We express the positions, kinetic energy, potential energy, and virtual power symbolically.
The time derivatives ``\dot{p}_c``, ``\dot{p}_p`` are computed automatically by `D.(...)`.

```@example main
p_c = [x, 0.0]
p_p = [x + l * sin(θ), l * cos(θ)]
F   = [u, 0.0]

T      = 0.5 * m_c * sum(D.(p_c) .^ 2) + 0.5 * m_p * sum(D.(p_p) .^ 2)
V      = g * (m_p * p_p[2])
W_ctrl = transpose(D.(p_c)) * F   # virtual power of the control force
```

### Euler–Lagrange Equations and Mass-Matrix Inversion

Starting from ``\mathcal{L} = T - V``, Symbolics.jl computes the three terms of the
Euler–Lagrange equations and assembles the residual vector. Substituting static aliases
``(v,\omega,\dot{v},\dot\omega)`` for the time derivatives makes it possible to identify
the mass matrix ``M`` as the Jacobian of the residual with respect to the accelerations
``(\dot{v}, \dot\omega)``. The system ``M\,a = -b`` is then solved analytically.

```@example main
L = T - V

A = D.(Symbolics.gradient(L, D.(q)))
B = Symbolics.gradient(L, q)
Q = Symbolics.gradient(W_ctrl, D.(q))   # generalised forces

# Euler-Lagrange residual: d/dt(∂L/∂q̇) - ∂L/∂q - Q = 0
el_eqs = expand_derivatives.(A - B - Q)

# Freeze time derivatives into static algebraic variables
sub_rules = Dict(D(x) => v, D(θ) => ω, D(D(x)) => dv, D(D(θ)) => dω)
res = Symbolics.substitute.(el_eqs, (sub_rules,))

# Identify mass matrix M and bias vector b such that M·[dv; dω] + b = 0
Mass = Symbolics.jacobian(res, [dv, dω])
bias = Symbolics.substitute.(res, (Dict(dv => 0.0, dω => 0.0),))

# Analytically invert the 2×2 mass matrix
accel = Symbolics.simplify_fractions.(Mass \ (-bias))

# Fully explicit state derivative: Ẋ = [v, ω, v̇, ω̇]
dx_dt = [v, ω, accel[1], accel[2]]
```

### Code Generation

`build_function` compiles the symbolic expression `dx_dt` into a native Julia function.
The `force_SA=true` flag generates a **StaticArrays** kernel, which avoids heap allocations
inside the ODE right-hand side — important for solver performance. Parameters are stored as
an `SVector` to match the generated kernel's expected input type.

```@example main
f_expr = build_function(dx_dt, [x, θ, v, ω], u, [m_c, m_p, l, g];
    expression=Val{false}, force_SA=true)
f_cartpole = f_expr[1]   # out-of-place variant: (state, u, params) → SVector

const p_vals = SA[m_c_val, m_p_val, l_val, g_val]   # SVector, matches SA kernel

cartpole_dynamics(X, U) = f_cartpole(X, U, p_vals)
```

### Optimal Control Problem Definition

We now formulate the optimal control problem using the `@def` macro from OptimalControl.jl.
The boundary conditions ``X(0) = X(t_f)`` encode the periodicity of the orbit directly.

```@example main
@def cartpole_ocp begin
    t ∈ [0, tf_val], time
    X = (x_ocp, θ_ocp, v_ocp, ω_ocp) ∈ R⁴, state
    F_ctrl ∈ R, control

    X(0)      == [0, 0, 0, 0.1]
    X(tf_val) == [0, 0, 0, 0.1]   # periodicity: X(tf) = X(0)

    Ẋ(t) == cartpole_dynamics(X(t), F_ctrl(t))

    ∫(F_ctrl(t)^2) → min
end
```

### Solving the NLP

The problem is transcribed into a nonlinear program using direct collocation on a uniform
grid of 100 intervals, then handed to **Ipopt** via NLPModelsIpopt. The `@time` macro
reports wall-clock time, which on first call includes Julia's JIT compilation.

```@example main
initial_guess = (state=[0.0, 0.0, 0.0, 0.1], control=0.0)

@time sol = solve(cartpole_ocp; display=true, grid_size=100, init=initial_guess)
```

### Results

```@example main
println("--- Optimal Limit Cycle Found ---")
println("Total Control Energy (∫F² dt): ", sol.objective)

tsol = time_grid(sol)
Xsol = state(sol).(tsol)
Fsol = control(sol).(tsol)

# Stack state vectors into a 4×N matrix and unpack rows
xsol, θsol, vsol, ωsol = eachrow(reduce(hcat, Xsol))

qplot  = plot(tsol, [xsol θsol], label=["x" "θ"],   title="Configuration")
dqplot = plot(tsol, [vsol ωsol], label=["v" "ω"],   title="Velocities")
uplot  = plot(tsol, Fsol,        label="u",          title="Control", linetype=:steppost)

plot(qplot, dqplot, uplot, layout=3, size=(800, 600))
```

The three panels show the cart position ``x`` and pendulum angle ``\theta``, the
corresponding velocities ``\dot{x}`` and ``\dot\theta``, and the optimal control force
``u``. The periodicity of the state trajectory confirms that a genuine limit cycle has been
found.

### Animation

The animation below shows the cart-pole evolving along the optimal limit-cycle trajectory.
The blue cart slides on the horizontal rail while the pendulum swings around the upright
equilibrium. A ghost trace follows the red bob to make the motion easier to read. Use the
controls to play, pause, or reset the animation and to adjust the playback speed.

```@example main
# Serialise trajectory for the JS animation (manual, no extra dependencies)
json_t  = "[" * join(string.(tsol),  ",") * "]"
json_x  = "[" * join(string.(xsol),  ",") * "]"
json_th = "[" * join(string.(θsol),  ",") * "]"

uid    = string(rand(UInt32), base=16)   # unique suffix → safe for multi-example pages
tf_str = string(round(tsol[end], digits=3))

html_str = """
<div style="text-align:center;font-family:sans-serif;margin:1.5em 0;">
  <canvas id="cpC$(uid)" width="700" height="380"
    style="border:1px solid #ccc;border-radius:6px;background:#fafafa;display:block;margin:0 auto;"></canvas>
  <div style="margin:6px 0;font-size:0.9em;color:#444;">
    <span id="cpT$(uid)">t = 0.000 / $(tf_str) s</span>
  </div>
  <div style="display:flex;justify-content:center;gap:10px;align-items:center;flex-wrap:wrap;margin-top:4px;">
    <button id="cpPl$(uid)" style="padding:5px 14px;cursor:pointer;">&#9654; Play</button>
    <button id="cpPa$(uid)" style="padding:5px 14px;cursor:pointer;">&#9646;&#9646; Pause</button>
    <button id="cpRe$(uid)" style="padding:5px 14px;cursor:pointer;">&#x21BA; Reset</button>
    <label style="font-size:0.9em;">Speed:
      <input type="range" id="cpSp$(uid)" min="0.25" max="3" step="0.25" value="1"
             style="vertical-align:middle;width:100px;">
      <span id="cpSv$(uid)">1.00&times;</span>
    </label>
  </div>
</div>
<script>
(function(){
  var T=$(json_t), X=$(json_x), TH=$(json_th), L=$(l_val);
  var cv=document.getElementById("cpC$(uid)");
  var ctx=cv.getContext("2d");
  var W=cv.width, H=cv.height;
  var railY=H*0.62, sc=80, cW=60, cH=28;
  var tSim=0, spd=1, rid=null, lts=null, ghost=[];
  function lerp(a,b,t){return a+(b-a)*t;}
  function interp(tq){
    var n=T.length;
    if(tq<=T[0])return{x:X[0],th:TH[0]};
    if(tq>=T[n-1])return{x:X[n-1],th:TH[n-1]};
    var lo=0,hi=n-2;
    while(lo<hi){var m=(lo+hi)>>1;if(T[m+1]<tq)lo=m+1;else hi=m;}
    var a=(tq-T[lo])/(T[lo+1]-T[lo]);
    return{x:lerp(X[lo],X[lo+1],a),th:lerp(TH[lo],TH[lo+1],a)};
  }
  function draw(tq){
    var d=interp(tq);
    var cx=W/2+d.x*sc, cy=railY;
    var px=cx, py=cy-cH/2;
    var bx=px+L*sc*Math.sin(d.th), by=py-L*sc*Math.cos(d.th);
    ghost.push({x:bx,y:by}); if(ghost.length>25)ghost.shift();
    ctx.clearRect(0,0,W,H);
    ctx.beginPath();ctx.moveTo(20,railY);ctx.lineTo(W-20,railY);
    ctx.strokeStyle="#aaa";ctx.lineWidth=3;ctx.stroke();
    var wr=8;
    [[cx-18,cy+cH/2+wr],[cx+18,cy+cH/2+wr]].forEach(function(p){
      ctx.beginPath();ctx.arc(p[0],p[1],wr,0,2*Math.PI);
      ctx.fillStyle="#888";ctx.fill();
    });
    ctx.beginPath();ctx.rect(cx-cW/2,cy-cH/2,cW,cH);
    ctx.fillStyle="#4a90d9";ctx.fill();
    ctx.strokeStyle="#2c5f8a";ctx.lineWidth=1.5;ctx.stroke();
    ghost.forEach(function(g,i){
      ctx.beginPath();ctx.arc(g.x,g.y,5,0,2*Math.PI);
      ctx.fillStyle="rgba(220,60,60,"+((i+1)/ghost.length*0.35)+")";ctx.fill();
    });
    ctx.beginPath();ctx.moveTo(px,py);ctx.lineTo(bx,by);
    ctx.strokeStyle="#222";ctx.lineWidth=3;ctx.stroke();
    ctx.beginPath();ctx.arc(px,py,5,0,2*Math.PI);ctx.fillStyle="#555";ctx.fill();
    ctx.beginPath();ctx.arc(bx,by,10,0,2*Math.PI);
    ctx.fillStyle="#dc3c3c";ctx.fill();ctx.strokeStyle="#8a1a1a";ctx.lineWidth=1.5;ctx.stroke();
    document.getElementById("cpT$(uid)").textContent=
      "t = "+tq.toFixed(3)+" / "+T[T.length-1].toFixed(2)+" s";
  }
  function step(ts){
    if(lts!==null){
      var dt=(ts-lts)*0.001*spd; tSim+=dt;
      var tf=T[T.length-1]; if(tSim>tf){tSim-=tf;ghost=[];}
    }
    lts=ts; draw(tSim); rid=requestAnimationFrame(step);
  }
  document.getElementById("cpPl$(uid)").onclick=function(){
    if(rid===null){lts=null;rid=requestAnimationFrame(step);}
  };
  document.getElementById("cpPa$(uid)").onclick=function(){
    if(rid!==null){cancelAnimationFrame(rid);rid=null;}
  };
  document.getElementById("cpRe$(uid)").onclick=function(){
    if(rid!==null){cancelAnimationFrame(rid);rid=null;}
    tSim=0;ghost=[];lts=null;draw(0);
  };
  document.getElementById("cpSp$(uid)").oninput=function(){
    spd=parseFloat(this.value);
    document.getElementById("cpSv$(uid)").textContent=spd.toFixed(2)+"\u00d7";
  };
  draw(0);
})();
</script>
"""

HTML(html_str)
```