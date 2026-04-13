# [Direct and indirect methods for the Goddard problem](@id tutorial-goddard)

```@meta
Draft = false
```

## Introduction

```@raw html
<img src="./assets/Goddard_and_Rocket.jpg" style="float: left; margin: auto 10px;" width="200px">
```

For this example, we consider the well-known Goddard problem[^1] [^2] which models the ascent of a rocket through the atmosphere, and we restrict here ourselves to vertical (one dimensional) trajectories. The state variables are the altitude $r$, speed $v$ and mass $m$ of the rocket during the flight, for a total dimension of 3. The rocket is subject to gravity $g$, thrust $u$ and drag force $D$ (function of speed and altitude). The final time $t_f$ is free, and the objective is to reach a maximal altitude with a bounded fuel consumption.

We thus want to solve the optimal control problem in Mayer form

```math
    r(t_f) \to \max
```

subject to the controlled dynamics

```math
    \dot{r} = v, \quad
    \dot{v} = \frac{T_{\max}\,u - D(r,v)}{m} - g, \quad
    \dot{m} = -u,
```

and subject to the control constraint $u(t) \in [0,1]$ and the state constraint $v(t) \leq v_{\max}$. The initial state is fixed while only the final mass is prescribed.

!!! note "Nota bene"

    The Hamiltonian is affine with respect to the control, so singular arcs may occur,
    as well as constrained arcs due to the path constraint on the velocity (see below).

We import the [OptimalControl.jl](https://control-toolbox.org/OptimalControl.jl) package to define the optimal control problem and [NLPModelsIpopt.jl](https://jso.dev/NLPModelsIpopt.jl) to solve it. We import the [Plots.jl](https://docs.juliaplots.org) package to plot the solution. The [OrdinaryDiffEq.jl](https://docs.sciml.ai/OrdinaryDiffEq) package is used to define the shooting function for the indirect method and the [MINPACK.jl](https://github.com/sglyon/MINPACK.jl) package permits to solve the shooting equation.

```@example main-goddard
using OptimalControl  # to define the optimal control problem and more
using NLPModelsIpopt  # to solve the problem via a direct method
using OrdinaryDiffEq  # to get the Flow function from OptimalControl
using MINPACK         # NLE solver: use to solve the shooting equation
using Plots           # to plot the solution
```

## Optimal control problem

We define the problem

```@example main-goddard
const t0 = 0      # initial time
const r0 = 1      # initial altitude
const v0 = 0      # initial speed
const m0 = 1      # initial mass
const vmax = 0.1  # maximal authorized speed
const mf = 0.6    # final mass to target

goddard = @def begin # definition of the optimal control problem

    tf ∈ R, variable
    t ∈ [t0, tf], time
    x = (r, v, m) ∈ R³, state
    u ∈ R, control

    x(t0) == [r0, v0, m0]
    m(tf) == mf,         (1)
    0 ≤ u(t) ≤ 1
    r(t) ≥ r0
    0 ≤ v(t) ≤ vmax

    ẋ(t) == F0(x(t)) + u(t) * F1(x(t))

    -r(tf) → min

end

# Dynamics
const Cd = 310
const Tmax = 3.5
const β = 500
const b = 2

F0(x) = begin
    r, v, m = x
    D = Cd * v^2 * exp(-β*(r - 1)) # Drag force
    return [v, -D/m - 1/r^2, 0]
end

F1(x) = begin
    r, v, m = x
    return [0, Tmax/m, -b*Tmax]
end
nothing # hide
```

## Direct method

We then solve it

```@example main-goddard
direct_sol = solve(goddard; grid_size=250)
nothing # hide
```

and plot the solution

```@example main-goddard
plt = plot(direct_sol, label="direct", size=(800, 800))
```

## [Structure of the solution](@id tutorial-goddard-structure)

We first determine visually the structure of the optimal solution which is composed of a bang arc with maximal control, followed by a singular arc, then by a boundary arc and the final arc is with zero control. In summary, the structure is: **bang** ($u=1$) → **singular** → **boundary** ($v=v_{\max}$) → **off** ($u=0$). Note that the switching function vanishes along the singular and boundary arcs.

```@example main-goddard
t = time_grid(direct_sol)   # the time grid as a vector
x = state(direct_sol)       # the state as a function of time
u = control(direct_sol)     # the control as a function of time
p = costate(direct_sol)     # the costate as a function of time

H1 = Lift(F1)           # H1(x, p) = p' * F1(x)
φ(t) = H1(x(t), p(t))   # switching function
g(x) = vmax - x[2]      # state constraint v ≤ vmax

u_plot  = plot(t, u,     linewidth=2, label = "u(t)")
H1_plot = plot(t, φ,     linewidth=2, label = "H₁(x(t), p(t))")
g_plot  = plot(t, g ∘ x, linewidth=2, label = "g(x(t))")

plot(u_plot, H1_plot, g_plot, layout=(3,1), size=(600, 600))
```

We are now in position to solve the problem by an indirect shooting method. For an introduction to the indirect simple shooting method, see the [Indirect simple shooting](@ref tutorial-indirect-simple-shooting) tutorial. We first define the four control laws in feedback form and their associated flows. For this we need to compute some Lie derivatives, namely [Poisson brackets](https://en.wikipedia.org/wiki/Poisson_bracket) of Hamiltonians (themselves obtained as lifts to the cotangent bundle of vector fields), or derivatives of functions along a vector field. For instance, the control along the *minimal order* singular arcs is obtained as the quotient

```math
u_s = -\frac{H_{001}}{H_{101}}
```

of length three Poisson brackets:

```math
H_{001} = \{H_0,\{H_0,H_1\}\}, \quad H_{101} = \{H_1,\{H_0,H_1\}\}
```

where, for two Hamiltonians $H$ and $G$,

```math
\{H,G\} := (\nabla_p H|\nabla_x G) - (\nabla_x H|\nabla_p G).
```

While the Lie derivative of a function $f$ *wrt.* a vector field $X$ is simply obtained as

```math
(X \cdot f)(x) := f'(x) \cdot X(x),
```

and is used to the compute the control along the boundary arc,

```math
u_b(x) = -(F_0 \cdot g)(x) / (F_1 \cdot g)(x),
```

as well as the associated multiplier for the *order one* state constraint on the velocity:

```math
\mu(x, p) = H_{01}(x, p) / (F_1 \cdot g)(x).
```

!!! note "Poisson bracket and Lie derivative"

    The Poisson bracket $\{H,G\}$ is also given by the Lie derivative of $G$ along the Hamiltonian vector field $X_H = (\nabla_p H, -\nabla_x H)$ of $H$, that is

    ```math
        \{H,G\} = X_H \cdot G
    ```

    which is the reason why we use the `@Lie` macro to compute Poisson brackets below.

With the help of differential geometry primitives, these expressions are straightforwardly translated into Julia code. For more details on the differential geometry tools, see [Differential geometry tools](@extref OptimalControl manual-differential-geometry). For a simpler example involving a minimal order singular arc, see [Singular control](@extref OptimalControl example-singular-control).

```@example main-goddard
# Controls
const u0 = 0                            # off control
const u1 = 1                            # bang control

H0 = Lift(F0)                           # H0(x, p) = p' * F0(x)
H01  = @Lie {H0, H1}
H001 = @Lie {H0, H01}
H101 = @Lie {H1, H01}
us(x, p) = -H001(x, p) / H101(x, p)     # singular control

ub(x) = -(F0⋅g)(x) / (F1⋅g)(x)          # boundary control
μ(x, p) = H01(x, p) / (F1⋅g)(x)         # multiplier associated to the state constraint g

# Flows
f0 = Flow(goddard, (x, p, tf) -> u0)
f1 = Flow(goddard, (x, p, tf) -> u1)
fs = Flow(goddard, (x, p, tf) -> us(x, p))
fb = Flow(goddard, (x, p, tf) -> ub(x), (x, u, tf) -> g(x), (x, p, tf) -> μ(x, p))
nothing # hide
```

## Shooting function

Then, we define the shooting function according to the optimal structure we have determined, that is a concatenation of four arcs. The shooting function has 7 unknowns (3 components of the initial costate `p0` and 4 times: `t1`, `t2`, `t3`, `tf`) and 7 equations.

```@example main-goddard
x0 = [r0, v0, m0] # initial state

function shoot!(s, p0, t1, t2, t3, tf)

    x1, p1 = f1(t0, x0, p0, t1)
    x2, p2 = fs(t1, x1, p1, t2)
    x3, p3 = fb(t2, x2, p2, t3)
    xf, pf = f0(t3, x3, p3, tf)

    s[1] = xf[3] - mf           # final mass constraint
    s[2:3] = pf[1:2] - [1, 0]   # transversality conditions: r(tf) and v(tf) are free
                                # objective is -r(tf) → p_r(tf) = 1, and v(tf) free → p_v(tf) = 0
    s[4] = H1(x1, p1)           # H1 = H01 = 0
    s[5] = H01(x1, p1)          # at the entrance of the singular arc
    s[6] = g(x2)                # g = 0 when entering the boundary arc
    s[7] = H0(xf, pf)           # since tf is free

end
nothing # hide
```

## Initial guess

To solve the problem by an indirect shooting method, we then need a good initial guess, that is a good approximation of the initial costate, the three switching times and the final time.

```@example main-goddard
η = 1e-3
t13 = t[ abs.(φ.(t)) .≤ η ]
t23 = t[ 0 .≤ (g ∘ x).(t) .≤ η ]
p0 = p(t0)
t1 = min(t13...)
t2 = min(t23...)
t3 = max(t23...)
tf = t[end]

println("p0 = ", p0)
println("t1 = ", t1)
println("t2 = ", t2)
println("t3 = ", t3)
println("tf = ", tf)

# Norm of the shooting function at solution
using LinearAlgebra: norm
s = similar(p0, 7)
shoot!(s, p0, t1, t2, t3, tf)
println("\nNorm of the shooting function: ‖s‖ = ", norm(s), "\n")
```

## Indirect shooting

We aggregate the data to define the initial guess vector.

```@example main-goddard
ξ_guess = [p0..., t1, t2, t3, tf] # initial guess
```

### MINPACK.jl

We can use the [MINPACK.jl](https://github.com/sglyon/MINPACK.jl) package to solve the shooting equation. To compute the Jacobian of the shooting function we use the [DifferentiationInterface.jl](https://juliadiff.org/DifferentiationInterface.jl/DifferentiationInterface) package with [ForwardDiff.jl](https://juliadiff.org/ForwardDiff.jl) backend.

```@setup main-goddard
using NonlinearSolve  # interface to NLE solvers
struct MYSOL
    x::Vector{Float64}
end
function fsolve(f, j, x; kwargs...)
    try
        MINPACK.fsolve(f, j, x; kwargs...)
    catch e
        println("Error using MINPACK")
        println(e)
        println("hybrj not supported. Replaced by NonlinearSolve even if it is not visible on the doc.")
        nle! = (s, ξ, _) -> f(s, ξ)
        prob = NonlinearProblem(nle!, x)
        sol = solve(prob; show_trace=Val(true))
        return MYSOL(sol.u)
    end
end
```

```@example main-goddard
using DifferentiationInterface
import ForwardDiff
backend = AutoForwardDiff()
nothing # hide
```

Let us define the problem to solve.

```@example main-goddard
# auxiliary function with aggregated inputs
shoot!(s, ξ) = shoot!(s, ξ[1:3], ξ[4], ξ[5], ξ[6], ξ[7])

# Jacobian of the (auxiliary) shooting function
jshoot!(js, ξ) = jacobian!(shoot!, similar(ξ), js, backend, ξ)
nothing # hide
```

We are now in position to solve the problem with the `hybrj` solver from MINPACK.jl through the `fsolve` function, providing the Jacobian. Let us solve the problem and retrieve the initial costate solution.

```@example main-goddard
# resolution of S(ξ) = 0
minpack_sol = fsolve(shoot!, jshoot!, ξ_guess, show_trace=true)

# we retrieve the costate solution together with the times
ξ = minpack_sol.x
p0, t1, t2, t3, tf = ξ[1:3], ξ[4:end]...

println("MINPACK results:")
println("p0 = ", p0)
println("t1 = ", t1)
println("t2 = ", t2)
println("t3 = ", t3)
println("tf = ", tf)

# Norm of the shooting function at solution
s = similar(p0, 7)
shoot!(s, p0, t1, t2, t3, tf)
println("\nNorm of the shooting function: ‖s‖ = ", norm(s), "\n")
@assert norm(s) < 1e-6 "Indirect shooting failed for MINPACK"
```

### NonlinearSolve.jl

Alternatively, we can use the [NonlinearSolve.jl](https://docs.sciml.ai/NonlinearSolve) package to solve the shooting equation. The code is similar, but we use the `solve` function instead of `fsolve`. Let us define the problem.

```@example main-goddard
shoot!(s, ξ, λ) = shoot!(s, ξ[1:3], ξ[4], ξ[5], ξ[6], ξ[7])
prob = NonlinearProblem(shoot!, ξ_guess)
nothing # hide
```

Now, let us solve the problem and retrieve the initial costate and times.

```@example main-goddard
# resolution of S(ξ) = 0
sciml_sol = solve(prob; show_trace=Val(true))

# we retrieve the costate solution together with the times
ξ = sciml_sol.u
p0, t1, t2, t3, tf = ξ[1:3], ξ[4:end]...

println("\nNonlinearSolve results:")
println("p0 = ", p0)
println("t1 = ", t1)
println("t2 = ", t2)
println("t3 = ", t3)
println("tf = ", tf)

# Norm of the shooting function at solution
s = similar(p0, 7)
shoot!(s, p0, t1, t2, t3, tf)
println("\nNorm of the shooting function: ‖s‖ = ", norm(s), "\n")
@assert norm(s) < 1e-6 "Indirect shooting failed for NonlinearSolve"
```

### Comparison of the two solvers

Let us compare the results obtained by the two solvers.

```@example main-goddard
ξ_minpack = minpack_sol.x
ξ_sciml = sciml_sol.u

println("Comparison of the two solvers:")
println("MINPACK solution: ξ = ", ξ_minpack)
println("SciML solution:   ξ = ", ξ_sciml)
println("Difference:       Δξ = ", ξ_minpack - ξ_sciml)
println("Relative error:   ‖Δξ‖/‖ξ‖ = ", norm(ξ_minpack - ξ_sciml) / norm(ξ_minpack))
```

### Benchmarking

The results found for by the two solvers are extremely close, so now, lets benchmark these two resolutions to compare their performances.

```@example main-goddard
using BenchmarkTools
```

```@example main-goddard
@benchmark fsolve(shoot!, jshoot!, ξ_guess; tol=1e-8, show_trace=false) #MINPACK
```

```@example main-goddard
@benchmark solve(prob; abstol=1e-8, reltol=1e-8, show_trace=Val(false)) # NonlinearSolve
```

According to the NonlinearSolve documentation, for small nonlinear systems, it could be faster to use the [`SimpleNewtonRaphson()` descent algorithm](https://docs.sciml.ai/NonlinearSolve/stable/tutorials/code_optimization/).

```@example main-goddard
@benchmark solve(prob, SimpleNewtonRaphson(); abstol=1e-8, reltol=1e-8, show_trace=Val(false)) # NonlinearSolve
```

There exist different alternatives to solve the shooting equation, as shown in the benchmarks above.

## [Plot of the solution](@id tutorial-goddard-plot)

We plot the solution of the indirect solution (in red) over the solution of the direct method (in blue).

```@example main-goddard
f = f1 * (t1, fs) * (t2, fb) * (t3, f0) # concatenation of the flows
indirect_sol = f((t0, tf), x0, p0)      # compute the solution: state, costate, control...

plot!(plt, indirect_sol; label="indirect", color=2)
```

!!! note "Numerical comparison between direct and indirect solutions"

    ```@raw html
    <details><summary>Click to unfold and get the code for numerical comparisons.</summary>
    ```

    ```@example main-goddard
    using Printf

    function L2_norm(T, X)
        # T and X are supposed to be one dimensional
        s = 0.0
        for i in 1:(length(T) - 1)
            s += 0.5 * (X[i]^2 + X[i + 1]^2) * (T[i + 1] - T[i])
        end
        return √(s)
    end

    function print_numerical_comparisons(direct_sol, indirect_sol)

        # get time grids
        t_dir = time_grid(direct_sol)
        t_ind = time_grid(indirect_sol)

        # create common time grid (union of both grids)
        t_common = unique(sort([t_dir..., t_ind...]))

        # interpolate both solutions on common grid
        x_dir_func = state(direct_sol)
        u_dir_func = control(direct_sol)
        x_ind_func = state(indirect_sol)
        u_ind_func = control(indirect_sol)

        x_dir = x_dir_func.(t_common)
        u_dir = u_dir_func.(t_common)
        x_ind = x_ind_func.(t_common)
        u_ind = u_ind_func.(t_common)

        v_dir = variable(direct_sol)
        o_dir = objective(direct_sol)
        i_dir = iterations(direct_sol)

        v_ind = variable(indirect_sol)
        o_ind = objective(indirect_sol)

        x_vars = ["r", "v", "m"]
        u_vars = ["u"]
        v_vars = ["tf"]

        println("┌─ Goddard problem: direct vs indirect")
        println("│")
        println("├─  Number of Iterations")
        @printf("│     Direct: %d\n", i_dir)
        println("│")

        # States
        println("├─  States (L2 Norms)")
        for i in eachindex(x_vars)
            xi_dir = [x_dir[k][i] for k in eachindex(t_common)]
            xi_ind = [x_ind[k][i] for k in eachindex(t_common)]
            L2_ae = L2_norm(t_common, xi_dir - xi_ind)
            L2_re = L2_ae / (0.5 * (L2_norm(t_common, xi_dir) + L2_norm(t_common, xi_ind)))
            @printf("│     %-6s Abs: %.3e   Rel: %.3e\n", x_vars[i], L2_ae, L2_re)
        end
        println("│")

        # Controls
        println("├─  Controls (L2 Norms)")
        for i in eachindex(u_vars)
            ui_dir = [u_dir[k][i] for k in eachindex(t_common)]
            ui_ind = [u_ind[k][i] for k in eachindex(t_common)]
            L2_ae = L2_norm(t_common, ui_dir - ui_ind)
            L2_re = L2_ae / (0.5 * (L2_norm(t_common, ui_dir) + L2_norm(t_common, ui_ind)))
            @printf("│     %-6s Abs: %.3e   Rel: %.3e\n", u_vars[i], L2_ae, L2_re)
        end
        println("│")

        # Variables
        println("├─  Variables")
        for i in eachindex(v_vars)
            vi_dir = v_dir[i]
            vi_ind = v_ind[i]
            vi_ae = abs(vi_dir - vi_ind)
            vi_re = vi_ae / (0.5 * (abs(vi_dir) + abs(vi_ind)))
            @printf("│     %-6s Abs: %.3e   Rel: %.3e\n", v_vars[i], vi_ae, vi_re)
        end
        println("│")

        # Objective
        println("├─  Objective")
        o_ae = abs(o_dir - o_ind)
        o_re = o_ae / (0.5 * (abs(o_dir) + abs(o_ind)))
        @printf("│            Abs: %.3e   Rel: %.3e\n", o_ae, o_re)
        println("└─")
        return nothing
    end
    nothing # hide
    ```

    ```@raw html
    </details>
    ```

We now compare numerically the direct and indirect solutions. The comparison includes L2 norms for the state and control variables, as well as absolute and relative errors for the variable (final time) and the objective. The L2 norm is computed using the trapezoidal rule over the time grid of the direct solution.

```@example main-goddard
print_numerical_comparisons(direct_sol, indirect_sol)
```

## Hamiltonian verification

The Hamiltonian should be constant along the optimal trajectory and equal to zero since the final time is free. Let us verify this for both the direct and indirect solutions.

```@example main-goddard
# Hamiltonian function
H(x, p, u) = H0(x, p) + u * H1(x, p)

# Direct solution
t_dir = time_grid(direct_sol)
x_dir = state(direct_sol)
u_dir = control(direct_sol)
p_dir = costate(direct_sol)
H_direct = [H(x_dir(ti), p_dir(ti), u_dir(ti)) for ti in t_dir]

# Indirect solution
t_ind = time_grid(indirect_sol)
x_ind = state(indirect_sol)
u_ind = control(indirect_sol)
p_ind = costate(indirect_sol)
H_indirect = [H(x_ind(ti), p_ind(ti), u_ind(ti)) for ti in t_ind]

println("Direct method:")
println("  H(t0) = ", H_direct[1])
println("  H variation: max|H(t) - H(t0)| = ", maximum(abs.(H_direct .- H_direct[1])))

println("\nIndirect method:")
println("  H(t0) = ", H_indirect[1])
println("  H variation: max|H(t) - H(t0)| = ", maximum(abs.(H_indirect .- H_indirect[1])))
```

We can also plot the Hamiltonian along the trajectory to verify its constancy visually.

```@example main-goddard
# Plot
plot(t_dir, H_direct, label="Direct", linewidth=2)
plot!(t_ind, H_indirect, label="Indirect", linewidth=2, linestyle=:dash)
xlabel!("Time")
ylabel!("H(t)")
title!("Hamiltonian along the trajectory")
```

## References

[^1]: R.H. Goddard. A Method of Reaching Extreme Altitudes, volume 71(2) of Smithsonian Miscellaneous Collections. Smithsonian institution, City of Washington, 1919.

[^2]: H. Seywald and E.M. Cliff. Goddard problem in presence of a dynamic pressure limit. Journal of Guidance, Control, and Dynamics, 16(4):776–781, 1993.
