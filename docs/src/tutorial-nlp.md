# [NLP and DOCP manipulations](@id tutorial-nlp)

```@meta
Draft = false
CurrentModule =  OptimalControl
```

We describe here some more advanced operations related to the discretized optimal control problem. When calling `solve(ocp)`, three steps are performed internally:

1. The OCP is discretized into a DOCP (discretized optimal control problem)
2. An NLP (nonlinear programming) model is built from the DOCP
3. The NLP is solved, and a functional solution of the OCP is rebuilt from the NLP solution

These steps can also be done manually, which gives you more control over the solving process. This is useful when you want to:

- Use your own NLP solver
- Access and inspect the NLP model
- Benchmark different components separately
- Fine-tune each step of the resolution

Let us load the packages.

```@example main-nlp
using OptimalControl
using Plots
```

We define a simple test problem: a double integrator with minimal control energy.

```@example main-nlp
ocp = @def begin

    t ∈ [0, 1], time
    x ∈ R², state
    u ∈ R, control

    x(0) == [ -1, 0 ]
    x(1) == [ 0, 0 ]

    ẋ(t) == [ x₂(t), u(t) ]

    ∫( 0.5u(t)^2 ) → min

end
nothing # hide
```

## Step-by-step resolution with OptimalControl

We now perform the resolution step by step, using OptimalControl strategy instances.

### Step 1: Build the initial guess

First, we build an initial guess for the problem. Here we pass `nothing` to use the default initialization.

```@example main-nlp
init = build_initial_guess(ocp, nothing)
nothing # hide
```

### Step 2: Discretize the OCP

We create a discretizer strategy and use it to discretize the optimal control problem into a DOCP.

```@example main-nlp
discretizer = OptimalControl.Collocation(grid_size=100, scheme=:trapeze)
docp = discretize(ocp, discretizer)
nothing # hide
```

The `docp` is a `DiscretizedModel` that contains information about the discretization, including a copy of the original OCP.

### Step 3: Build the NLP model

Next, we create a modeler strategy and build the NLP model from the discretized problem.

```@example main-nlp
modeler = OptimalControl.ADNLP(backend=:optimized)
nlp = nlp_model(docp, init, modeler)
nothing # hide
```

The `nlp` is an `ADNLPModel` (from the NLPModels ecosystem) representing the discretized nonlinear programming problem.

### Step 4: Solve the NLP

We have two approaches to solve the NLP problem.

#### Approach A: Using OptimalControl solver wrapper

We can create an OptimalControl solver strategy and use it to solve the NLP:

```@example main-nlp
solver = OptimalControl.Ipopt(print_level=5, tol=1e-8, mu_strategy="adaptive")
nlp_sol = solve(nlp, solver; display=true)
nothing # hide
```

#### Approach B: Using external NLP solvers directly

Alternatively, we can use NLP solvers directly from their respective packages. For instance, with [NLPModelsIpopt.jl](https://jso.dev/NLPModelsIpopt.jl):

```@example main-nlp
using NLPModelsIpopt
nlp_sol_ipopt = ipopt(nlp; print_level=5, mu_strategy="adaptive", tol=1e-8, sb="yes")
nothing # hide
```

Or with [MadNLP.jl](https://madnlp.github.io/MadNLP.jl):

```@example main-nlp
using MadNLP
nlp_sol_madnlp = madnlp(nlp; print_level=MadNLP.ERROR, tol=1e-8)
nothing # hide
```

### Step 5: Build the OCP solution

Finally, we build the optimal control solution from the NLP solution and plot it. Note that the multipliers from the NLP solver are used to compute the costate.

```@example main-nlp
sol = ocp_solution(docp, nlp_sol, modeler)
plot(sol)
```

## All-in-one with explicit strategies

All the steps above can also be performed in a single `solve` call by passing the strategy instances explicitly:

```@example main-nlp
sol = solve(ocp, init, discretizer, modeler, solver; display=false)
plot(sol)
```

This is the explicit mode of `solve`, which gives you full control over each component while keeping the interface simple. See the [explicit mode documentation](@extref OptimalControl manual-solve-explicit) for more details.

## Initial guess and warm start

An initial guess can be passed at different stages of the process.

If you have a previous solution, you can use it to warm-start a new resolution:

```@example main-nlp
# Use the previous solution as initial guess
init_warm = build_initial_guess(ocp, sol)

# Discretize and solve with warm start
docp_warm = discretize(ocp, discretizer)
nlp_warm = nlp_model(docp_warm, init_warm, modeler)
nlp_sol_warm = solve(nlp_warm, solver; display=false)
sol_warm = ocp_solution(docp_warm, nlp_sol_warm, modeler)

println("Iterations with warm start: ", iterations(sol_warm))
println("Iterations without warm start: ", iterations(sol))
```

You can also pass the initial guess directly to `solve` in explicit mode:

```@example main-nlp
sol_warm2 = solve(ocp, init_warm, discretizer, modeler, solver; display=false)
println("Objective value: ", objective(sol_warm2))
```
