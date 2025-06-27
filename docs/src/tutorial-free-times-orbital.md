```@meta
Draft = false
```

# [Optimal control problem with free final orbital time](@id tutorial-free-times-orbital)




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
```
And plot the solution.

```@example orbit
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