# [NLP and DOCP manipulations](@id tutorial-nlp)

```@meta
CurrentModule =  OptimalControl
```

We describe here some more advanced operations related to the discretized optimal control problem.
When calling `solve(ocp)` three steps are performed internally:

- first, the OCP is discretized into a DOCP (a nonlinear optimization problem),
- then, this DOCP is solved with a nonlinear programming (NLP) solver, which returns a solution of the discretized problem,
- finally, a functional solution of the OCP is rebuilt from the solution of the discretized problem.

These steps can also be done separately, for instance if you want to use your own NLP solver. 

Let us load the packages.

```@example main-nlp
using OptimalControl
using Plots
```

We define a test problem

```@example main-nlp
ocp = @def begin

    t ∈ [0, 1], time
    x ∈ R², state
    u ∈ R, control

    x(0) == [ -1, 0 ]
    x(1) == [ 0, 0 ]

    ẋ(t) == [ x₂(t), u(t) ]

    ∫( 0.5u(t)^2 ) → min

end
nothing # hide
```

## Discretization and NLP problem

We discretize the problem.

```@example main-nlp
docp, nlp = direct_transcription(ocp)
nothing # hide
```

The DOCP contains information related to the transcription, including a copy of the original OCP, and the NLP is the resulting discretized problem, in our case an `ADNLPModel`.

We can now use the solver of our choice to solve it.

## Resolution of the NLP problem

For a first example we use the `ipopt` solver from [NLPModelsIpopt.jl](https://jso.dev/NLPModelsIpopt.jl) package to solve the NLP problem.

```@example main-nlp
using NLPModelsIpopt
nlp_sol = ipopt(nlp; print_level=5, mu_strategy="adaptive", tol=1e-8, sb="yes")
nothing # hide
```

Then we can rebuild and plot an optimal control problem solution (note that the multipliers are optional, but the OCP costate will not be retrieved if the multipliers are not provided).

```@example main-nlp
sol = build_OCP_solution(docp; primal=nlp_sol.solution, dual=nlp_sol.multipliers)
plot(sol)
```
## Change the NLP solver

Alternatively, we can use [MadNLP.jl](https://madnlp.github.io/MadNLP.jl) to solve anew the NLP problem:

```@example main-nlp
using MadNLP
nlp_sol = madnlp(nlp)
```

## Initial guess

An initial guess, including warm start, can be passed to [`direct_transcription`](https://control-toolbox.org/OptimalControl.jl/stable/dev-ctdirect.html#CTDirect.direct_transcription-Tuple{Model,%20Vararg{Any}}) the same way as for `solve`.

```@example main-nlp
docp, nlp = direct_transcription(ocp; init=sol)
nothing # hide
```

It can also be changed after the transcription is done, with  [`set_initial_guess`](https://control-toolbox.org/OptimalControl.jl/stable/dev-ctdirect.html#CTDirect.set_initial_guess-Tuple{CTDirect.DOCP,%20Any,%20Any}).

```@example main-nlp
set_initial_guess(docp, nlp, sol)
nothing # hide
```
