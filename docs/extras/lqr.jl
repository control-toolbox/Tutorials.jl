using Pkg
Pkg.activate("docs")

using OptimalControl
using NLPModelsIpopt
using Plots

function lqr(tf)
    x0 = [
        0
        1
    ]

    A = [
        0 1
        -1 0
    ]

    B = [
        0
        1
    ]

    Q = [
        1 0
        0 1
    ]

    R = 1

    ocp = @def begin
        t ∈ [0, tf], time
        x ∈ R², state
        u ∈ R, control
        x(0) == x0
        ẋ(t) == [x₂(t), - x₁(t) + u(t)]
        0.5∫(x₁(t)^2 + x₂(t)^2 + u(t)^2) → min
    end

    return ocp
end

solutions = []   # empty list of solutions
tfs = [3, 5, 30]

for tf in tfs
    solution = solve(lqr(tf); display=false)
    push!(solutions, solution)
end

plt = plot(solutions[1]; time=:normalize)
for sol in solutions[2:end]
    plot!(plt, sol; time=:normalize)
end

# we plot only the state and control variables and we add the legend
N = length(tfs)
px1 = plot(plt[1]; legend=false, xlabel="s", ylabel="x₁")
px2 = plot(
    plt[2]; label=reshape(["tf = $tf" for tf in tfs], (1, N)), xlabel="s", ylabel="x₂"
)
pu = plot(plt[5]; legend=false, xlabel="s", ylabel="u")

using Plots.PlotMeasures # for leftmargin, bottommargin
plot(px1, px2, pu; layout=(1, 3), size=(800, 300), leftmargin=5mm, bottommargin=5mm)

####

plt = plot(solutions[1], :state, :control; time=:normalize, label="tf = $(tfs[1])")
for (tf, sol) in zip(tfs[2:end], solutions[2:end])
    plot!(plt, sol, :state, :control; time=:normalize, label="tf = $tf")
end

px1 = plot(plt[1]; legend=false, xlabel="s", ylabel="x₁")
px2 = plot(plt[2]; legend=true, xlabel="s", ylabel="x₂")
pu = plot(plt[3]; legend=false, xlabel="s", ylabel="u")
plot(px1, px2, pu; layout=(1, 3), size=(800, 300), leftmargin=5mm, bottommargin=5mm)
