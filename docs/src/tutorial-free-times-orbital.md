# [Optimal control problem with free final orbital time](@id tutorial-free-times-orbital)

```@meta
Draft = true
```

## A more concrete example about the change of orbit of a satellite:

## Here are the required packages for the tutorial:

```@example orbit
using OptimalControl
using NLPModelsIpopt
using Plots
using Printf
using LinearAlgebra
using NLsolve
```

## Definition of the problem

```@example orbit

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
        ẋ(t) == [
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

function optimal_control(x, p, tf)
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
# Create flow using OptimalControl.jl
flow = Flow(ocp, optimal_control)
```


## Shooting function

The shooting function encodes the boundary conditions:

```@example orbit
function shooting_function!(s, ξ)
    p0 = ξ[1:4]
    tf = ξ[5]

    x0 = [x₁₀, x₂₀, x₃₀, x₄₀]

    try
        # Use OptimalControl.jl Flow interface with variable
        x_final, p_final = flow(0.0, x0, p0, tf, tf)

        x₁f, x₂f, x₃f, x₄f = x_final
        p₁f, p₂f, p₃f, p₄f = p_final

        s[1] = x₁f^2 + x₂f^2 - r_f^2
        s[2] = p₁f
        s[3] = p₂f

        u_final = optimal_control(x_final, p_final, tf)
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
    
    # Define time points for plotting FIRST
    t_indirect = range(0, tf_opt, length=1000)
    
    # Solve the optimal trajectory
    sol_indirect = flow((0.0, tf_opt), x0, p0_opt)
    
    # Extract trajectories by evaluating at time points
    x_traj = zeros(length(t_indirect), 4)
    p_traj = zeros(length(t_indirect), 4)
    u_traj = zeros(length(t_indirect), 2)
    
    for (i, t) in enumerate(t_indirect)
        x_traj[i, :], p_traj[i, :] = flow(0.0, x0, p0_opt, t)
        u_traj[i, :] = optimal_control(x_traj[i, :], p_traj[i, :], tf_opt)    
    end

    println("Indirect solution generated successfully!")
    println("Final time: $(tf_opt)")
    println("Final position: $(x_traj[end, 1:2])")
    println("Final radius: $(sqrt(x_traj[end, 1]^2 + x_traj[end, 2]^2))")
end
```
# Visualisation of results
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
    
    t_indirect = range(0, tf_opt, length=1000)
    
    x_traj = zeros(length(t_indirect), 4)
    p_traj = zeros(length(t_indirect), 4)
    u_traj = zeros(length(t_indirect), 2)
    
    # Use flow function directly to evaluate at each time point
    for (i, t) in enumerate(t_indirect)
        x_traj[i, :], p_traj[i, :] = flow(0.0, x0, p0_opt, t)
        u_traj[i, :] = optimal_control(x_traj[i, :], p_traj[i, :], tf_opt)
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