# Navigation problem, MPC approach

```@meta
Draft = false
```

We consider a ship in a constant current $w=(w_x,w_y)$, where $\|w\|<1$. 
The [heading angle](https://en.wikipedia.org/wiki/Heading) is controlled, leading to the following differential equations:

```math
\begin{array}{rcl}
    \dot{x}(t) &=& w_x + \cos\theta(t), \quad t \in [0, t_f] \\
    \dot{y}(t) &=& w_y + \sin\theta(t), \\
    \dot{\theta}(t) &=& u(t). 
\end{array}
```

The angular velocity is limited and normalized: $\|u(t)\| \leq 1$. There are boundary conditions at the initial time $t=0$ and at the final time $t=t_f$, on the position $(x,y)$ and on the angle $\theta$. The objective is to minimize the final time. This topic stems from a collaboration between the University Côte d'Azur and the French company [CGG](https://www.cgg.com), which is interested in optimal maneuvers of very large ships for marine exploration.

```@raw html
<img 
      src="./assets/ship.jpg"
      style="display: block;
      margin-left: auto;
      margin-right: auto;
      width: 50%;"
>
```

## Data 

```@example main
using OptimalControl, NLPModelsIpopt, Plots, OrdinaryDiffEq, LinearAlgebra, Plots.PlotMeasures

t0 = 0.
x0 = 0. 
y0 = 0.
θ0 = π/7
xf = 4.
yf = 7.
θf = -π/2

function current(x, y) # current as a function of position
    ε = 1e-1
    w = [ 0.6, 0.4 ]
    δw = ε * [ y, -x ] / sqrt(1+x^2+y^2)
    w = w + δw
    if (w[1]^2 + w[2]^2 >= 1)
        error("|w| >= 1")
    end
    return w
end

#
function plot_state!(plt, x, y, θ; color=1)
    plot!(plt, [x], [y], marker=:circle, legend=false, color=color, markerstrokecolor=color, markersize=5, z_order=:front)
    quiver!(plt, [x], [y], quiver=([cos(θ)], [sin(θ)]), color=color, linewidth=2, z_order=:front)
    return plt
end

function plot_current!(plt; current=current, N=10, scaling=1)
    for x ∈ range(xlims(plt)..., N)
        for y ∈ range(ylims(plt)..., N)
            w = scaling*current(x, y)
            quiver!(plt, [x], [y], quiver=([w[1]], [w[2]]), color=:black, linewidth=0.5, z_order=:back)
        end
    end
    return plt
end

# Display the boundary conditions and the current in the augmented phase plane
plt = plot(
    xlims=(-2, 6), 
    ylims=(-1, 8), 
    size=(600, 600), 
    aspect_ratio=1, 
    xlabel="x", 
    ylabel="y", 
    title="Boundary Conditions",
    leftmargin=5mm, 
    bottommargin=5mm,
)

plot_state!(plt, x0, y0, θ0; color=2)
plot_state!(plt, xf, yf, θf; color=2)
annotate!([(x0, y0, ("q₀", 12, :top)), (xf, yf, ("qf", 12, :bottom))])
plot_current!(plt)
```

```@example main
function plot_trajectory!(plt, t, x, y, θ; N=5) # N: number of points where we will display θ

    # trajectory
    plot!(plt, x.(t), y.(t), legend=false, color=1, linewidth=2, z_order=:front)

    if N > 0

        # length of the path
        s = 0
        for i ∈ 2:length(t)
            s += norm([x(t[i]), y(t[i])] - [x(t[i-1]), y(t[i-1])])
        end

        # interval of length
        Δs = s/(N+1)
        tis = []
        s = 0
        for i ∈ 2:length(t)
            s += norm([x(t[i]), y(t[i])] - [x(t[i-1]), y(t[i-1])])
            if s > Δs && length(tis) < N
                push!(tis, t[i])
                s = 0
            end
        end

        # display intermediate points
        for ti ∈ tis
            plot_state!(plt, x(ti), y(ti), θ(ti); color=1)
        end

    end

    return plt
    
end
nothing # hide
```

## OptimalControl solver

```@example main
function solve(t0, x0, y0, θ0, xf, yf, θf, w; 
    grid_size=300, tol=1e-8, max_iter=500, print_level=4, display=true)

    # Definition of the problem
    ocp = @def begin

        tf ∈ R, variable
        t ∈ [t0, tf], time
        q = (x, y, θ) ∈ R³, state
        u ∈ R, control

        -1 ≤ u(t) ≤ 1

        -2 ≤ x(t) ≤ 6
        -2 ≤ y(t) ≤ 8
        -2π ≤ x(t) ≤ 2π

        q(t0) == [x0, y0, θ0]
        q(tf) == [xf, yf, θf]

        q̇(t) == [w[1]+cos(θ(t)), 
                  w[2]+sin(θ(t)), 
                  u(t)]

        tf → min

    end

    # Initialization
    tf_init = 1.5*norm([xf, yf]-[x0, y0])
    x_init(t) = [ x0, y0, θ0 ] * (tf_init-t)/(tf_init-t0) + [xf, yf, θf] * (t-t0)/(tf_init-t0)
    u_init = (θf - θ0) / (tf_init-t0)
    init = (state=x_init, control=u_init, variable=tf_init)

    # Resolution
    sol = OptimalControl.solve(ocp; 
        init=init,
        grid_size=grid_size, 
        tol=tol, 
        max_iter=max_iter, 
        print_level=print_level, 
        display=display, 
        disc_method=:euler,
    )

    # Retrieval of useful data
    t = time_grid(sol)
    q = state(sol)
    x = t -> q(t)[1]
    y = t -> q(t)[2]
    θ = t -> q(t)[3]
    u = control(sol)
    tf = variable(sol)
    
    return t, x, y, θ, u, tf, iterations(sol), sol.solver_infos.constraints_violation
    
end
nothing # hide
```

## First resolution

We consider a constant current and we solve a first time the problem.

```@example main
# Resolution
t, x, y, θ, u, tf, iter, cons = solve(t0, x0, y0, θ0, xf, yf, θf, current(x0, y0); display=false);

println("Iterations: ", iter)
println("Constraints violation: ", cons)
println("tf: ", tf)

# Displaying the trajectory
plt_q = plot(xlims=(-2, 6), ylims=(-1, 8), aspect_ratio=1, xlabel="x", ylabel="y")
plot_state!(plt_q, x0, y0, θ0; color=2)
plot_state!(plt_q, xf, yf, θf; color=2)
plot_current!(plt_q; current=(x, y) -> current(x0, y0))
plot_trajectory!(plt_q, t, x, y, θ)

# Displaying the control
plt_u = plot(t, u; color=1, legend=false, linewidth=2, xlabel="t", ylabel="u")

# Final display
plot(plt_q, plt_u; 
    layout=(1, 2), 
    size=(1200, 600),
    leftmargin=5mm, 
    bottommargin=5mm,
    plot_title="Constant Current Simulation"
)
```

## Simulation of the Real System

In the previous simulation, we assumed that the current is constant. However, from a practical standpoint, the current depends on the position $(x, y)$. Given a current model, provided by the function `current`, we can simulate the actual trajectory of the ship, as long as we have the initial condition and the control over time.

```@example main
function realistic_trajectory(tf, t0, x0, y0, θ0, u, current; abstol=1e-12, reltol=1e-12, saveat=[])
    
    function rhs!(dq, q, dummy, t)
        x, y, θ = q
        w = current(x, y)
        dq[1] = w[1] + cos(θ)
        dq[2] = w[2] + sin(θ)
        dq[3] = u(t)
    end
    
    q0 = [x0, y0, θ0]
    tspan = (t0, tf)
    ode = ODEProblem(rhs!, q0, tspan)
    sol = OrdinaryDiffEq.solve(ode, Tsit5(), abstol=abstol, reltol=reltol, saveat=saveat)

    t = sol.t
    x = t -> sol(t)[1]
    y = t -> sol(t)[2]
    θ = t -> sol(t)[3]

    return t, x, y, θ
    
end
nothing # hide
```

```@example main
# Realistic trajectory
t, x, y, θ = realistic_trajectory(tf, t0, x0, y0, θ0, u, current)

# Displaying the trajectory
plt_q = plot(xlims=(-2, 6), ylims=(-1, 8), aspect_ratio=1, xlabel="x", ylabel="y")
plot_state!(plt_q, x0, y0, θ0; color=2)
plot_state!(plt_q, xf, yf, θf; color=2)
plot_current!(plt_q; current=current)
plot_trajectory!(plt_q, t, x, y, θ)
plot_state!(plt_q, x(tf), y(tf), θ(tf); color=3)

# Displaying the control
plt_u = plot(t, u; color=1, legend=false, linewidth=2, xlabel="t", ylabel="u")

# Final display
plot(plt_q, plt_u; 
    layout=(1, 2), 
    size=(1200, 600),
    leftmargin=5mm, 
    bottommargin=5mm,
    plot_title="Simulation with Current Model"
)
```

## MPC Approach

In practice, we do not have the actual current data for the entire trajectory in advance, which is why we will regularly recalculate the optimal control. The idea is to update the optimal control at regular time intervals, taking into account the current at the position where the ship is located. We are therefore led to solve a number of problems with constant current, with this being updated regularly. This is an introduction to the so-called Model Predictive Control (MPC) methods.

```@example main
function MPC(t0, x0, y0, θ0, xf, yf, θf, current)

    Nmax = 20   # maximum number of iterations for the MPC method
    ε = 1e-1    # radius on the final condition to stop calculations
    Δt = 1.0    # fixed time step for the MPC method
    P = 300      # number of discretization points for the solver

    t1 = t0
    x1 = x0
    y1 = y0
    θ1 = θ0

    data = []

    N = 1
    stop = false

    while !stop
        
        # Retrieve the current at the current position
        w = current(x1, y1)

        # Solve the problem
        t, x, y, θ, u, tf, iter, cons = solve(t1, x1, y1, θ1, xf, yf, θf, w; grid_size=P, display=false);

        # Calculate the next time
        if (t1 + Δt < tf)
            t2 = t1 + Δt
        else
            t2 = tf
            println("t2=tf: ", t2)
            stop = true
        end

        # Store the data: the current initial time, the next time, the control
        push!(data, (t2, t1, x(t1), y(t1), θ(t1), u, tf))

        # Update the parameters of the MPC method: simulate reality
        t, x, y, θ = realistic_trajectory(t2, t1, x1, y1, θ1, u, current)
        t1 = t2
        x1 = x(t1)
        y1 = y(t1)
        θ1 = θ(t1)

        # Calculate the distance to the target position
        distance = norm([x1, y1, θ1] - [xf, yf, θf])
        println("N: ", N, "\t distance: ", distance, "\t iterations: ", iter, "\t constraints: ", cons, "\t tf: ", tf)
        if !((distance > ε) && (N < Nmax))
            stop = true
        end

        #
        N += 1

    end

    return data
end

data = MPC(t0, x0, y0, θ0, xf, yf, θf, current)
nothing # hide
```

## Display

```@example main
# Trajectory
plt_q = plot(xlims=(-2, 6), ylims=(-1, 8), aspect_ratio=1, xlabel="x", ylabel="y")

# Final condition
plot_state!(plt_q, xf, yf, θf; color=2)

# Current
plot_current!(plt_q; current=current)

# Control
plt_u = plot(xlabel="t", ylabel="u")

for d ∈ data

    t2, t1, x1, y1, θ1, u, tf = d

    # Calculate the actual trajectory
    t, x, y, θ = realistic_trajectory(t2, t1, x1, y1, θ1, u, current)

    # Trajectory
    plot_state!(plt_q, x1, y1, θ1; color=2)
    plot_trajectory!(plt_q, t, x, y, θ; N=0)

    # Control
    plot!(plt_u, t, u; color=1, legend=false, linewidth=2)

end

# last point
d = data[end]
t2, t1, x1, y1, θ1, u, tf = d
t, x, y, θ = realistic_trajectory(t2, t1, x1, y1, θ1, u, current)
plot_state!(plt_q, x(tf), y(tf), θ(tf); color=3)

#
plot(plt_q, plt_u; 
    layout=(1, 2), 
    size=(1200, 600),
    leftmargin=5mm, 
    bottommargin=5mm,
    plot_title="Simulation with Current Model"
)
```