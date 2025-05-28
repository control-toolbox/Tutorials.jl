using OptimalControl    # to define the optimal control problem and its flow
using OrdinaryDiffEq    # to get the Flow function from OptimalControl
using MINPACK           # NLE solver: use to solve the shooting equation
using Plots             # to plot the solution

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

u(x, p) = p
φ = Flow(ocp, u)

π((x, p)) = x

S(p0) = π( φ(t0, x0, p0, tf) ) - xf    # shooting function

ξ = [0.1]    # initial guess

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

using DifferentiationInterface
import ForwardDiff
backend = AutoForwardDiff()

nle!  = ( s, ξ) -> s[1] = S(ξ[1])                                 # auxiliary function
jnle! = (js, ξ) -> jacobian!(nle!, similar(ξ), js, backend, ξ)    # 


indirect_sol = fsolve(nle!, jnle!, ξ; show_trace=true)    # resolution of S(p0) = 0
p0_sol = indirect_sol.x[1]                                # costate solution
println("\ncostate:    p0 = ", p0_sol)
println("shoot: |S(p0)| = ", abs(S(p0_sol)), "\n")

sol = φ((t0, tf), x0, p0_sol)
plot(sol)

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

p0_sol

plot(S, p0_sol; size=(800, 450))