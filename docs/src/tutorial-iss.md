# [Indirect simple shooting](@id tutorial-indirect-simple-shooting)

In this tutorial we present the indirect simple shooting method on a simple example.

Let us start by importing the necessary packages. We import the [OptimalControl.jl](https://control-toolbox.org/OptimalControl.jl) package to define the optimal control problem.  We import the [Plots.jl](https://docs.juliaplots.org) package to plot the solution.  The [OrdinaryDiffEq.jl](https://docs.sciml.ai/OrdinaryDiffEq) package is used to define the shooting function for the indirect method and the [MINPACK.jl](https://github.com/sglyon/MINPACK.jl) package permits to solve the shooting equation.


```@example main-iss
using BenchmarkTools    # to benchmark the methods
using OptimalControl    # to define the optimal control problem and its flow
using OrdinaryDiffEq    # to get the Flow function from OptimalControl
using MINPACK           # NLE solver: use to solve the shooting equation
using NonlinearSolve    # interface to NLE solvers
using Plots             # to plot the solution
```

## Optimal control problem

Let us consider the following optimal control problem:

```math
\left\{ 
    \begin{array}{l}
        \min \displaystyle \frac{1}{2} \int_{t_0}^{t_f} u^2(t) \, \mathrm{d} t\\[1.0em]
        \dot{x}(t)  =  \displaystyle -x(t)+\alpha x^2(t)+u(t), \quad  u(t) \in \R, 
        \quad t \in [t_0, t_f] \text{ a.e.}, \\[0.5em]
        x(t_0) = x_0, \quad x(t_f) = x_f,
    \end{array}
\right.%
```

with $t_0 = 0$, $t_f = 1$, $x_0 = -1$, $x_f = 0$, $\alpha=1.5$ and $\forall\, t \in [t_0, t_f]$, $x(t) \in \R$.

```@example main-iss
const t0 = 0
const tf = 1
const x0 = -1
const xf = 0
const α  = 1.5
ocp = @def begin

    t ∈ [t0, tf], time
    x ∈ R, state
    u ∈ R, control

    x(t0) == x0
    x(tf) == xf

    ẋ(t) == -x(t) + α * x(t)^2 + u(t)

    ∫( 0.5u(t)^2 ) → min
    
end
nothing # hide
```

## Boundary value problem

The **pseudo-Hamiltonian** of this problem is

```math
    H(x,p,u) = p \, (-x+\alpha x^2+u) + p^0 u^2 /2,
```

where $p^0 = -1$ since we are in the normal case. From the Pontryagin Maximum Principle, the maximising control is given by

```math
u(x, p) = p
```

since $\partial^2_{uu} H = p^0 = - 1 < 0$. Plugging this control in feedback form into the pseudo-Hamiltonian, and considering the limit conditions, we obtain the following two-points boundary value problem (BVP).

```math
    \left\{ 
        \begin{array}{l}
            \dot{x}(t)  = \phantom{-} \nabla_p H[t] = -x(t) + \alpha x^2(t) + u(x(t), p(t)) 
            = -x(t) + \alpha x^2(t) + p(t), \\[0.5em]
            \dot{p}(t)  = -           \nabla_x H[t] = (1 - 2 \alpha x(t))\, p(t),    \\[0.5em]
            x(t_0)        = x_0, \quad x(t_f) = x_f,
        \end{array}
    \right.
```

where $[t]~=  (x(t),p(t),u(x(t), p(t)))$.

!!! note "Our goal"

    Our goal is to solve this boundary value problem (BVP), which is equivalent to solving the Pontryagin Maximum Principle (PMP), which provides necessary conditions for optimality.

## Shooting function

To achive our goal, let us first introduce the pseudo-Hamiltonian vector field

```math
    \vec{H}(z,u) = \left( \nabla_p H(z,u), -\nabla_x H(z,u) \right), \quad z = (x,p),
```

and then denote by $\varphi_{t_0, x_0, p_0}(\cdot)$ the solution of the following Cauchy problem

```math
\dot{z}(t) = \vec{H}(z(t), u(z(t))), \quad z(t_0) = (x_0, p_0).
```

Our goal becomes to solve

```math
\pi( \varphi_{t_0, x_0, p_0}(t_f) ) = x_f,
```

where $\pi(x, p) = x$. To compute $\varphi$ with [OptimalControl.jl](https://control-toolbox.org/OptimalControl.jl) package, we define the flow of the associated Hamiltonian vector field by:

```@example main-iss
u(x, p) = p
φ = Flow(ocp, u)
nothing # hide
```

We define also the projection function on the state space.

```@example main-iss
π((x, p)) = x
nothing # hide
```

!!! note "Nota bene"

    Actually, $\varphi_{t_0, x_0, p_0}(\cdot)$ is also solution of
    
    ```math
        \dot{z}(t) = \vec{\mathbf{H}}(z(t)), \quad z(t_0) = (x_0, p_0),
    ```
    where $\mathbf{H}(z) = H(z, u(z))$ and $\vec{\mathbf{H}} = (\nabla_p \mathbf{H}, -\nabla_x \mathbf{H})$. This is what is actually computed by `Flow`.

Now, to solve the (BVP) we introduce the **shooting function**:

```math
    \begin{array}{rlll}
        S \colon    & \R    & \longrightarrow   & \R \\
                    & p_0    & \longmapsto     & S(p_0) = \pi( \varphi_{t_0, x_0, p_0}(t_f) ) - x_f.
    \end{array}
```

```@example main-iss
S(p0) = π( φ(t0, x0, p0, tf) ) - xf    # shooting function
nothing # hide
```

## Resolution of the shooting equation

At the end, solving (BVP) is equivalent to solve $S(p_0) = 0$. This is what we call the **indirect simple shooting method**. We define an initial guess.

```@example main-iss
ξ = [0.1]    # initial guess
nothing # hide
```

### MINPACK.jl

We can use the [MINPACK.jl](https://github.com/sglyon/MINPACK.jl) to solve the shooting equation. To compute the Jacobian of the shooting function we use [DifferentiationInterface.jl](https://juliadiff.org/DifferentiationInterface.jl/DifferentiationInterface) with [ForwardDiff.jl](https://juliadiff.org/ForwardDiff.jl) backend.

```@setup main-iss
using MINPACK
function fsolve(f, j, x; kwargs...)
    try
        MINPACK.fsolve(f, j, x; kwargs...)
    catch e
        println("Erreur using MINPACK")
        println(e)
        println("hybrj not supported. Replaced by hybrd even if it is not visible on the doc.")
        MINPACK.fsolve(f, x; kwargs...)
    end
end
```

```@example main-iss
using DifferentiationInterface
import ForwardDiff
backend = AutoForwardDiff()
nothing # hide
```

Let us define the problem to solve.

```@example main-iss
nle!(s, ξ) = s[1] = S(ξ[1])                                 # auxiliary function
jnle!(js, ξ) = jacobian!(nle!, similar(ξ), js, backend, ξ)  # Jacobian of nle
nothing # hide
```

We are now in position to solve the problem with the `hybrj` solver from MINPACK.jl through the `fsolve` function, providing the Jacobian. Let us solve the problem and retrieve the initial costate solution.

```@example main-iss
sol = fsolve(nle!, jnle!, ξ; show_trace=true)    # resolution of S(p0) = 0
p0_sol = sol.x[1]                                # costate solution
println("\ncostate:    p0 = ", p0_sol)
println("shoot: |S(p0)| = ", abs(S(p0_sol)), "\n")
```

### NonlinearSolve.jl

Alternatively, we can use the [NonlinearSolve.jl](https://docs.sciml.ai/NonlinearSolve) package to solve the shooting equation. The code is similar, but we use the `solve` function instead of `fsolve`. Let us define the problem.

```@example main-iss
nle!(s, ξ, λ) = s[1] = S(ξ[1])    # auxiliary function
prob = NonlinearProblem(nle!, ξ)  # NLE problem with initial guess
nothing # hide
```

Now, let us solve the problem and retrieve the initial costate solution.

```@example main-iss
sol = solve(prob; show_trace=Val(true)) # resolution of S(p0) = 0  
p0_sol = sol.u[1] # costate solution
println("\ncostate:    p0 = ", p0_sol)
println("shoot: |S(p0)| = ", abs(S(p0_sol)), "\n")
```

### Benchmarking

Let us benchmark the methods to solve the shooting equation.

```@example main-iss
@benchmark fsolve(nle!, jnle!, ξ; tol=1e-8, show_trace=false) # MINPACK
```

```@example main-iss
@benchmark solve(prob; abstol=1e-8, reltol=1e-8, show_trace=Val(false)) # NonlinearSolve
```

According to the NonlinearSolve documentation, for small nonlinear systems, it could be faster to use the [`SimpleNewtonRaphson()` descent algorithm](https://docs.sciml.ai/NonlinearSolve/stable/tutorials/code_optimization/). 

```@example main-iss
@benchmark solve(prob, SimpleNewtonRaphson(); abstol=1e-8, reltol=1e-8, show_trace=Val(false)) # NonlinearSolve
```

## Plot of the solution

The solution can be plot calling first the flow.

```@example main-iss
sol = φ((t0, tf), x0, p0_sol)
plot(sol)
```

In the indirect shooting method, the search for the optimal control is replaced by the computation of its associated extremal. This computation is equivalent to finding the initial costate (or covector) that solves the shooting function. Let us now plot the extremal trajectory in the phase space, along with the shooting function and its solution.

```@raw html
<details><summary>Code of the plot function in the phase space.</summary>
```

```@example main-iss
using Plots.PlotMeasures 
function Plots.plot(S::Function, p0::Float64; Np0=20, kwargs...) 
 
    # times for wavefronts
    times = range(t0, tf, length=3)

    # times for trajectories
    tspan = range(t0, tf, length=100)

    # interval of initial covector
    p0_min = -0.5 
    p0_max = 2 

    # covector solution
    p0_sol = p0 
 
    # plot of the flow in phase space
    plt_flow = plot() 
    p0s = range(p0_min, p0_max, length=Np0) 
    for i ∈ eachindex(p0s) 
        sol = φ((t0, tf), x0, p0s[i])
        x = state(sol).(tspan)
        p = costate(sol).(tspan)
        label = i==1 ? "extremals" : false 
        plot!(plt_flow, x, p, color=:blue, label=label) 
    end 
 
    # plot of wavefronts in phase space 
    p0s = range(p0_min, p0_max, length=200) 
    xs  = zeros(length(p0s), length(times)) 
    ps  = zeros(length(p0s), length(times)) 
    for i ∈ eachindex(p0s) 
        sol = φ((t0, tf), x0, p0s[i], saveat=times)
        xs[i, :] .= state(sol).(times) 
        ps[i, :] .= costate(sol).(times) 
    end 
    for j ∈ eachindex(times) 
        label = j==1 ? "flow at times" : false 
        plot!(plt_flow, xs[:, j], ps[:, j], color=:green, linewidth=2, label=label) 
    end 
 
    #  
    plot!(plt_flow, xlims=(-1.1, 1), ylims=(p0_min, p0_max)) 
    plot!(plt_flow, [0, 0], [p0_min, p0_max], color=:black, xlabel="x", ylabel="p", label="x=xf") 
     
    # solution 
    sol = φ((t0, tf), x0, p0_sol)
    x = state(sol).(tspan)
    p = costate(sol).(tspan)
    plot!(plt_flow, x, p, color=:red, linewidth=2, label="extremal solution") 
    plot!(plt_flow, [x[end]], [p[end]], seriestype=:scatter, color=:green, label=false) 
 
    # plot of the shooting function  
    p0s = range(p0_min, p0_max, length=200) 
    plt_shoot = plot(xlims=(p0_min, p0_max), ylims=(-2, 4), xlabel="p₀", ylabel="y") 
    plot!(plt_shoot, p0s, S, linewidth=2, label="S(p₀)", color=:green) 
    plot!(plt_shoot, [p0_min, p0_max], [0, 0], color=:black, label="y=0") 
    plot!(plt_shoot, [p0_sol, p0_sol], [-2, 0], color=:black, label="p₀ solution", linestyle=:dash) 
    plot!(plt_shoot, [p0_sol], [0], seriestype=:scatter, color=:green, label=false) 
 
    # final plot 
    plot(plt_flow, plt_shoot; layout=(1,2), leftmargin=15px, bottommargin=15px, kwargs...) 
 
end
nothing # hide
```

```@raw html
</details>
```

```@example main-iss
plot(S, p0_sol; size=(800, 450))
```