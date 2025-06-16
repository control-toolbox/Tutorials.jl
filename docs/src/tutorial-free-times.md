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

# Verification of results
```@raw html
<details style="margin-left:3em"><summary>Verification of results.</summary>
```

Here the theorical part :
```math
H = p1*x2 + p2*u
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
```@raw html
</details>
```


Now we can compare the results found with the direct method whith the theoritical analysis :
```@example main-disc
tf = variable(sol)
u = control(sol)
p = costate(sol)
x = state(sol)
p° = -1

xf = x(tf)
pf = p(tf)

Htf = pf[1]*xf[2] + pf[2]*u(tf)
@printf("H(tf) = %.3f\n", Htf)
@printf("x(tf) = [%.3f, %.3f]\n", xf[1], xf[2])
@printf("p(tf) = [%.3f, %.3f]\n", pf[1], pf[2])
```

The numerical results closely match the theoretical predictions: the final state x(tf)=[1,0] is exactly satisfied.
The costate and Hamiltonian values at final time show a small deviation (≈ 0.01), likely due to numerical precision.
Overall, the direct method confirms the theoretical analysis with excellent accuracy.

We can analyse the influence of using different discretization sizes (grid_size), and observed the following results for the optimal tf:

```@example main-disc
for N in [20, 50, 100, 200]
    solN = solve(ocp; grid_size=N, display=false)
    @printf("grid_size = %3d → tf = %.5f\n", N, objective(solN))
end
```

This example shows that problems with a free final time can be sensitive to discretization. A small grid may lead to suboptimal or slightly inaccurate results.

## Indirect method

We keep the structure of the solution found with the direct method
```@example main-disc
t = time_grid(sol)
x = state(sol)
u = control(sol)
p = costate(sol)

H(x,p,u,t) = p[1](t)*x[2](t) + p[2](t)*u(u)
H(xf, pf, tf) = pf[1]*xf[2] + pf[2]*u(tf)

# Hamiltonian vector field
F1(x) = [0, 1]
H1 = Lift(F1)
φ(t) = H1(x(t), p(t))  # switching function
```


First, lets define the different flows
```@example main-disc
using OrdinaryDiffEq  # to get the Flow function from OptimalControl

const u1 = 1
const u0 = -1

f1 = Flow(ocp, (x, p, tf) -> u1)
f0 = Flow(ocp, (x, p, tf) -> u0)
```

And with this we have the following shooting function
```@example main-disc
x0 = [0.0, 0.0]  # état initial
xf_target = [1.0, 0.0]  # état final

function shoot!(s, p0, t1, tf)
x1, p1 = f1(0.0, x0, p0, t1)
xf, pf = f0(t1, x1, p1, tf)

# Conditions de tir
s[1:2] = xf .- xf_target  # état final
s[3]   = H(xf, pf, tf) - 1  # condition de transversalité H(tf) = 1
s[4]   = H1(x1, p1)       # φ = 0 au switch
end
```

Before solving our problem we must find a good initial guess to help the convergence of the algorithm
```@example main-disc
p0 = p(0.0)
φ_arr = abs.(φ.(t))
t1 = t[argmin(φ_arr)]
tf = t[end]

println("p0 ≈ ", p0)
println("t1 ≈ ", t1)
println("tf ≈ ", tf)

s = zeros(4)
shoot!(s, p0, t1, tf)

ξ = [p0..., t1, tf]
```

And we finally can solve the problem using an indirect method
```@example main-disc
using NonlinearSolve
using DifferentiationInterface
import ForwardDiff

backend = AutoForwardDiff()

struct MYSOL
    x::Vector{Float64}
end

function fsolve(f, j, x; kwargs...)
    try
        MINPACK.fsolve(f, j, x; kwargs...)
    catch e
        println("MINPACK error:")
        println(e)
        println("→ Using NonlinearSolve fallback")

        # Wrap pour respecter l'interface de NonlinearProblem
        function wrapped_f(s, ξ, p)
            f(s, ξ)
            return nothing
        end

        prob = NonlinearProblem(wrapped_f, x)
        sol = solve(prob; abstol=1e-8, reltol=1e-8)
        return MYSOL(sol.u)
    end
end

# Agrégation du problème
nle!  = (s, ξ) -> shoot!(s, ξ[1:2], ξ[3], ξ[4])
jnle! = (js, ξ) -> jacobian!(nle!, similar(ξ), js, backend, ξ)

ξ0 = [p0... , t1, tf]
indirect_sol = fsolve(nle!, jnle!, ξ0)

p0 = indirect_sol.x[1:2]
t1 = indirect_sol.x[3]
tf = indirect_sol.x[4]

f = f1 * (t1, f0)
flow_sol = f((0.0, tf), x0, p0)

plot!(flow_sol, label="indirect", color=:red)
```

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



## Verification of results
```@raw html
<details style="margin-left:3em"><summary>Verification of results.</summary>
```
Here is the theoretical part using Pontryagin's Maximum Principle:
```math
H = p_1 x_2 + p_2 u + 1
```

Conditions from Pontryagin’s theorem:
```math
p_1' = 0 \quad \Rightarrow \quad p_1 = c_1 \quad (\text{constant}) \\
p_2' = -p_1 \quad \Rightarrow \quad p_2 = -c_1 t + c_2
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

On \( t \in [t_0, t_s] \) :
```math
x_2' = u = 1 \quad \Rightarrow \quad x_2(t) = t - t_0 \\
x_1' = x_2 \quad \Rightarrow \quad x_1(t) = \frac{(t - t_0)^2}{2}
```

At switching time \( t = t_s \) :
```math
x_2(t_s) = t_s - t_0 \\
x_1(t_s) = \frac{(t_s - t_0)^2}{2}
```

On \( t \in [t_s, 0] \) :
```math
x_2' = u = -1 \quad \Rightarrow \quad x_2(t) = x_2(t_s) - (t - t_s) \\
x_1' = x_2 \quad \Rightarrow \quad x_1(t) = x_1(t_s) + \int_{t_s}^t x_2(s) ds
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
- Switching time: \( t_s = -1 \)
- Initial time: \( t_0 = -2 \)

Control:
```math
u(t) = 1 \quad \text{on} \quad [-2, -1] \\
u(t) = -1 \quad \text{on} \quad [-1, 0]
```
```@raw html
</details>
```
## An example with free final and free initial time :
```@example both_time
using OptimalControl
using NLPModelsIpopt
using BenchmarkTools
using DataFrames
using Plots
using Printf

function double_integrator_freet0tf()
    @def ocp begin
        v=(t0, tf) ∈ R², variable
        t ∈ [t0, tf], time
        x ∈ R², state
        u ∈ R, control
        -1 ≤ u(t) ≤ 1
        x(t0) == [0, 0]
        x(tf) == [1, 0]
        0.05 ≤ t0 ≤ 10
        0.05 ≤ tf ≤ 10
        0.01 ≤ tf - t0 ≤ Inf
        ẋ(t) == [x₂(t), u(t)]
        t0 → max
    end

    return ocp
end
nothing # hide
```

#  Direct resolution with both free times :

We now solve the problem using a direct method, with automatic treatment of the free initial time.

```@example both_time
ocp = double_integrator_freet0tf()
sol = solve(ocp; grid_size=100)
plot(sol; label="direct", size=(800, 800))
```


## A more concrete example about the change of orbit of a satellite :

```@example orbit
using OptimalControl
using NLPModelsIpopt
using BenchmarkTools
using DataFrames
using Plots
using Printf

const x₁₀ = -42272.67       # initial position x
const x₂₀ = 0               # initial position y
const x₃₀ = 0               # initial velocity in x
const x₄₀ = -5696.72        # initial velocity in y
const μ = 5.1658620912*1e12 # gravitational parameter
const γ_max = 0.05          # maximal thrust norm
const r_f = 42165             # target orbit radius (final distance to origin)


function min_orbit_tf()
    @def ocp begin
        tf ∈ R, variable
        t ∈ [0, tf], time
        x ∈ R⁴, state
        u ∈ R², control
        u₁(t)^2 + u₂(t)^2 ≤ γ_max^2
        x(0) == [x₁₀, x₂₀, x₃₀, x₄₀]
        x₁(tf)^2 + x₂(tf)^2 == r_f^2
        0.05 ≤ tf ≤ Inf
        ẋ(t) == [
            x₃(t),
            x₄(t),
            -μ * x₁(t) / ((x₁(t)^2 + x₂(t)^2)^(3/2)) + u₁(t),
            -μ * x₂(t) / ((x₁(t)^2 + x₂(t)^2)^(3/2)) + u₂(t)
        ]
        tf → min
    end

    return ocp
end
nothing # hide
```

#  Direct resolution with minimal orbital time :

We now solve the problem using a direct method, with automatic treatment of the free initial time.

```@example orbit
ocp = min_orbit_tf()
sol = solve(ocp;init=(variable=13.4,), grid_size=100)
plot(sol; label="direct", size=(800, 800))
```
