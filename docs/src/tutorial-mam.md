# Minimal Action Method using Optimal Control

```@meta
Draft = false
```

The Minimal Action Method (MAM) is a numerical technique for finding the most probable transition pathway between stable states in stochastic dynamical systems. It achieves this by minimizing an action functional that represents the path's deviation from the deterministic dynamics, effectively identifying the path of least resistance through the system's landscape.

This tutorial demonstrates how to implement MAM as an optimal control problem, using the classical Maier-Stein model as a benchmark example.

## Required Packages

```@example main-mam
using OptimalControl
using NLPModelsIpopt
using Plots, Printf
```

## Problem Statement

We aim to find the most probable transition path between two stable states of a stochastic dynamical system. For a system with deterministic dynamics $f(x)$ and small noise, the transition path minimizes the action functional:

```math
S[x(\cdot), u(\cdot)] = \int_0^T \|u(t) - f(x(t))\|^2 \, dt
```

subject to the path dynamics:

```math
\dot{x}(t) = u(t), \quad x(0) = x_0, \quad x(T) = x_f
```

where $x_0$ and $x_f$ are the initial and final states, and $T$ is the transition time.

!!! note "Physical interpretation"

    The action $S$ measures the "cost" of deviating from the deterministic flow $f(x)$. Paths with smaller action are exponentially more likely in the small noise limit.

## Problem Setup

We consider a 2D system with a double-well flow, called the Maier-Stein model. It is a famous benchmark problem as it exhibits non-gradient dynamics with two stable equilibrium points at $(-1,0)$ and $(1,0)$, connected by a non-trivial transition path.
The system's deterministic dynamics are given by:

```@example main-mam
# Define the vector field
f(u, v) = [u - u^3 - 10*u*v^2,  -(1 - u^2)*v]
f(x) = f(x...)
nothing # hide
```

## Optimal Control Formulation

The minimal action path minimizes the deviation from the deterministic dynamics:

```@example main-mam
function ocp(T)
    action = @def begin
        t ∈ [0, T], time
        x ∈ R², state
        u ∈ R², control
        x(0) == [-1, 0]                      # Starting point (left well)
        x(T) == [1, 0]                       # End point (right well)
        ẋ(t) == u(t)                         # Path dynamics
        ∫( sum((u(t) - f(x(t))).^2) ) → min  # Minimize deviation from deterministic flow
    end
    return action
end
nothing # hide
```

## Initial Guess

We provide an initial guess for the path using a simple interpolation with the `@init` macro:

```@example main-mam
# Time horizon
T = 50

# Helper functions for initial state guess
L(t) = -(1 - t/T) + t/T      # Linear interpolation from -1 to 1
P(t) = 0.3*(-L(t)^2 + 1)     # Parabolic arc (approximates saddle crossing)

init = @init ocp(T) begin
    # Linear interpolation for x₁
    x₁(t) := L(t)
    # Parabolic guess for x₂
    x₂(t) := P(t)
    # Control from vector field
    u(t) := f(L(t), P(t))
end
nothing # hide
```

!!! note "Initial guess strategy"

    The initial guess uses a simple geometric path: linear interpolation in $x_1$ and a parabolic arc in $x_2$. This provides a reasonable starting point that avoids the unstable saddle point at the origin. The control is initialized to follow the deterministic flow along this path.

## Solving the Problem

We solve the problem in two steps for better accuracy:

!!! note "Two-step resolution"

    Starting with a coarse grid (50 points) allows for faster initial convergence. Refining with a fine grid (1000 points) then improves accuracy of the solution.

```@example main-mam
# First solve with coarse grid
sol = solve(ocp(T); init=init, grid_size=50)

# Refine solution with finer grid
sol = solve(ocp(T); init=sol, grid_size=1000)

# Objective value
objective(sol)
```

## Visualizing Results

Let's plot the solution trajectory and phase space:

```@example main-mam
plot(sol)
```

```@example main-mam
# Phase space plot
MLP = state(sol).(time_grid(sol))
scatter(first.(MLP), last.(MLP), 
        title="Minimal Action Path",
        xlabel="u",
        ylabel="v",
        label="Transition path")
```

The resulting path shows the most likely transition between the two stable states given a transient time $T=50$, minimizing the action functional while respecting the system's dynamics.

## Minimize with respect to T

To find the maximum likelihood path, we also need to minimize the transient time `T`. Hence, we perform a discrete continuation over the parameter `T` by solving the optimal control problem over a continuous range of final times `T`, using each solution to initialize the next problem.

```@example main-mam
# Continuation function to avoid global variables
function continuation_mam(Ts; init_guess=init)
    objectives = Float64[]
    iterations_list = Int[]
    current_sol = init_guess
    
    for T in Ts
        current_sol = solve(ocp(T); display=false, init=current_sol, grid_size=1000, tol=1e-8)
        push!(objectives, objective(current_sol))
        push!(iterations_list, iterations(current_sol))
    end
    
    return objectives, iterations_list, current_sol
end

# Perform continuation
Ts = range(1, 100, 100)
objectives, iters, final_sol = continuation_mam(Ts)

# Display results
println(" Time   Objective     Iterations")
for i in eachindex(Ts)
    @printf("%6.2f  %9.6e  %d\n", Ts[i], objectives[i], iters[i])
end
nothing # hide
```

We can now analyze the results and identify the optimal transition time:

```@example main-mam
# Find optimal time
idx_min = argmin(objectives)
T_min = Ts[idx_min]
obj_min = objectives[idx_min]

@printf("Optimal transition time: T* = %.2f\n", T_min)
@printf("Minimal action: S* = %.6e\n", obj_min)
```

Let us visualize the evolution of the objective function with respect to the transition time:

```@example main-mam
plt1 = scatter(Ts, log10.(objectives), xlabel="Time", label="Objective (log10)")
vline!(plt1, [T_min], label="Minimum at T=$(round(T_min, digits=1))", z_order=:back)
plot(plt1, size=(800,400))
```

We now focus on the region around the minimum for a clearer view:

```@example main-mam
plt2 = scatter(Ts[40:100], log10.(objectives[40:100]), xlabel="Time", label="Objective (log10)")
vline!(plt2, [T_min], label="Minimum", z_order=:back)
plot(plt2, size=(800,400))
```

!!! note "Interpretation"

    The optimal transition time $T^*$ balances two competing effects: shorter times require larger deviations from the deterministic flow (higher action), while longer times allow the system to follow the flow more closely. The minimum represents the most probable transition time in the small noise limit.
