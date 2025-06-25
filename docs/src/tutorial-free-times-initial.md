```@meta
Draft = true
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
<details>
<summary> Click to show/hide mathematical verification.</summary>
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

#  Direct resolution with free initial time :

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
<details>
<summary> Click to show/hide mathematical verification.</summary>
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
## An example with both final and inital times being free:
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

#  Direct resolution with both final and inital times being free:


We now solve the problem using a direct method, with automatic treatment of the free initial time.

```@example both_time
ocp = double_integrator_freet0tf()
sol = solve(ocp; grid_size=100)
plot(sol; label="direct", size=(800, 800))
```


## A more concrete example about the change of orbit of a satellite:

```@example orbit
using OptimalControl
using NLPModelsIpopt
using BenchmarkTools
using DataFrames
using Plots
using Printf
using OrdinaryDiffEq
using NonlinearSolve
using DifferentiationInterface
import ForwardDiff
using LinearAlgebra
using NLsolve


const x₁₀ = -42272.67       # initial position x
const x₂₀ = 0               # initial position y
const x₃₀ = 0               # initial velocity in x
const x₄₀ = -5696.72        # initial velocity in y
const μ = 5.1658620912*1e12 # gravitational parameter
const γ_max = 0.05          # maximal thrust norm
const r_f = 42165             # target orbit radius (final distance to origin)
const rf3    = r_f^3
const α  = sqrt(μ/rf3);

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
# Indirect resolution with minimal orbital time :


## Hamiltonian and optimal control

From the Pontryagin Maximum Principle, the optimal control is bang-bang:

```@example orbit
ocp = min_orbit_tf()

function optimal_control(p)
    p₃, p₄ = p[3], p[4]
    p_norm = sqrt(p₃^2 + p₄^2)
    
    if p_norm > 1e-10
        u₁ = γ_max * p₃ / p_norm
        u₂ = γ_max * p₄ / p_norm
        return [u₁, u₂]
    else
        return [0.0, 0.0]
    end
end

function hamiltonian(x, p, u)
    x₁, x₂, x₃, x₄ = x
    p₁, p₂, p₃, p₄ = p
    u₁, u₂ = u
    
    r = sqrt(x₁^2 + x₂^2)
    r³ = r^3
    
    return (p₁*x₃ + p₂*x₄ + 
            p₃*(-μ*x₁/r³ + u₁) + 
            p₄*(-μ*x₂/r³ + u₂))
end
```

## Augmented dynamics
```@raw html
<details>
<summary>Click to show/hide Augmented dynamics code</summary>
```
We define the combined state-costate system:

```@example orbit
function augmented_dynamics!(dx, x_aug, params, t)
    x = x_aug[1:4]
    p = x_aug[5:8]
    
    x₁, x₂, x₃, x₄ = x
    p₁, p₂, p₃, p₄ = p
    
    r = sqrt(x₁^2 + x₂^2)
    r³ = r^3
    r⁵ = r^5
    
    u = optimal_control(p)
    u₁, u₂ = u
    
    # State dynamics (unchanged)
    dx[1] = x₃
    dx[2] = x₄
    dx[3] = -μ*x₁/r³ + u₁
    dx[4] = -μ*x₂/r³ + u₂
    
    # CORRECTED Costate dynamics: ṗ = -∂H/∂x
    dx[5] = -(p₃ * μ * (3*x₁^2/r⁵ - 1/r³) + p₄ * μ * (3*x₁*x₂/r⁵))  # -∂H/∂x₁
    dx[6] = -(p₃ * μ * (3*x₁*x₂/r⁵) + p₄ * μ * (3*x₂^2/r⁵ - 1/r³))   # -∂H/∂x₂
    dx[7] = -p₁  # -∂H/∂x₃
    dx[8] = -p₂  # -∂H/∂x₄
end

function create_flow()
    function flow_func(t_span, x0, p0)
        x_aug0 = vcat(x0, p0)
        prob = ODEProblem(augmented_dynamics!, x_aug0, t_span)
        sol = solve(prob, Tsit5(), reltol=1e-12, abstol=1e-12)
        return sol
    end
    return flow_func
end

flow = create_flow()
```
```@raw html
</details>
```
## Shooting function

The shooting function encodes the boundary conditions:

```@example orbit
function shooting_function!(s, ξ)
    p0 = ξ[1:4]
    tf = ξ[5]

    x0 = [x₁₀, x₂₀, x₃₀, x₄₀]

    try
        sol = flow((0.0, tf), x0, p0)

        x_final = sol.u[end][1:4]
        p_final = sol.u[end][5:8]

        x₁f, x₂f, x₃f, x₄f = x_final
        p₁f, p₂f, p₃f, p₄f = p_final

        s[1] = x₁f^2 + x₂f^2 - r_f^2
        s[2] = p₁f
        s[3] = p₂f

        u_final = optimal_control(p_final)
        H_final = hamiltonian(x_final, p_final, u_final)
        s[4] = H_final

        s[5] = p₁f^2 + p₂f^2 + p₃f^2 + p₄f^2 - 1.0

    catch e
        fill!(s, 1e6)
    end

    return nothing
end
```

## Initial guess and solver setup

```@example orbit
# Initial guess for the shooting variables [p₀, tf]

function get_initial_guess()
    tf_guess = 13.4
    
    p0_guess = [1.0323e-4, 4.915e-5, 3.568e-4, -1.554e-4]
    
    return vcat(p0_guess, tf_guess)
end

ξ_init = get_initial_guess()
println("Initial guess: p₀ = $(ξ_init[1:4]), tf = $(ξ_init[5])")
```

## Solve the shooting problem

```@example orbit
# Solve the shooting problem using NLsolve
function solve_indirect_method()
    println("Solving indirect method...")

    # Set up the shooting problem
    result = nlsolve(shooting_function!, ξ_init,
                    method=:trust_region,
                    ftol=1e-10,
                    xtol=1e-10,
                    iterations=1000,
                    show_trace=false)

    if result.f_converged || result.x_converged
        println("Shooting method converged!")
        println("Final residual norm: $(norm(result.zero))")
        return result.zero
    else
        println("Shooting method failed to converge")
        println("Final residual norm: $(norm(result.zero))")
        return nothing
    end
end

ξ_solution = solve_indirect_method()

if ξ_solution !== nothing
    p0_opt = ξ_solution[1:4]
    tf_opt = ξ_solution[5]

    println("\nOptimal solution:")
    println("p₀ = $(p0_opt)")
    println("tf = $(tf_opt)")

    # Verify the solution
    s_check = zeros(5)
    shooting_function!(s_check, ξ_solution)
    println("Shooting function residual: $(s_check)")
end
```

## Generate the indirect solution trajectory

```@example orbit
if ξ_solution !== nothing
    # Generate the optimal trajectory
    x0 = [x₁₀, x₂₀, x₃₀, x₄₀]
    p0_opt = ξ_solution[1:4]
    tf_opt = ξ_solution[5]
    
    # Solve the optimal trajectory
    sol_indirect = flow((0.0, tf_opt), x0, p0_opt)
    
    # Extract time points for plotting
    t_indirect = range(0, tf_opt, length=1000)
    
    # Evaluate solution at time points
    x_traj = zeros(length(t_indirect), 4)
    p_traj = zeros(length(t_indirect), 4)
    u_traj = zeros(length(t_indirect), 2)
    
    for (i, t) in enumerate(t_indirect)
        state_full = sol_indirect(t)
        x_traj[i, :] = state_full[1:4]
        p_traj[i, :] = state_full[5:8]
        u_traj[i, :] = optimal_control(state_full[5:8])
    end
    
    println("Indirect solution generated successfully!")
    println("Final time: $(tf_opt)")
    println("Final position: $(x_traj[end, 1:2])")
    println("Final radius: $(sqrt(x_traj[end, 1]^2 + x_traj[end, 2]^2))")
end
```
# Visualistion of results
```@raw html
<details>
<summary>Click to show/hide indirect method visualization code</summary>
```

```@example orbit


# Simple visualization - just the basic plots
if ξ_solution !== nothing
    p0_opt = ξ_solution[1:4]
    tf_opt = ξ_solution[5]
    
    x0 = [x₁₀, x₂₀, x₃₀, x₄₀]
    sol_indirect = flow((0.0, tf_opt), x0, p0_opt)
    
    t_indirect = range(0, tf_opt, length=1000)
    
    x_traj = zeros(length(t_indirect), 4)
    p_traj = zeros(length(t_indirect), 4)
    u_traj = zeros(length(t_indirect), 2)
    
    for (i, t) in enumerate(t_indirect)
        state_full = sol_indirect(t)
        x_traj[i, :] = state_full[1:4]
        p_traj[i, :] = state_full[5:8]
        u_traj[i, :] = optimal_control(state_full[5:8])
    end
    
    # Combined plot with 3 rows
    plt = plot(layout=(3,1), size=(800, 900))
    
    plot!(plt[1], t_indirect, [x_traj[:,1] x_traj[:,2] x_traj[:,3] x_traj[:,4]], 
          label=["x₁" "x₂" "x₃" "x₄"], title="States")
    
    plot!(plt[2], t_indirect, [p_traj[:,1] p_traj[:,2] p_traj[:,3] p_traj[:,4]], 
          label=["p₁" "p₂" "p₃" "p₄"], title="Costates")
    
    plot!(plt[3], t_indirect, [u_traj[:,1] u_traj[:,2]], 
          label=["u₁" "u₂"], title="Control")
    
    plt  # This returns the plot for documentation
end
```
```@raw html
</details>
```