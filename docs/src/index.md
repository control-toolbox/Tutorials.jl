# Getting Started

This collection of tutorials is part of the [control-toolbox ecosystem](https://github.com/control-toolbox). The control-toolbox ecosystem gathers Julia packages for mathematical control and applications. It aims to provide tools to model and solve optimal control problems with ordinary differential equations by direct and indirect methods, both on CPU and GPU. If you want to define an optimal control problem and solve it, please check the [documentation](https://control-toolbox.org/OptimalControl.jl).

From this page, you can find a list of tutorials to solve optimal control problems with OptimalControl.

## Available tutorials

| Tutorial | Description |
| :------- | :---------- |
| **[Discrete continuation](@ref tutorial-continuation)** | Uses warm start to implement a discrete continuation method, solving a sequence of OCPs where each solution initializes the next. |
| **[Discretisation methods](@ref tutorial-discretisation-methods)** | Compares discretisation schemes (trapeze, midpoint, Gauss-Legendre…) and automatic differentiation backends on the Goddard problem. |
| **[Free final time](@ref tutorial-free-times-final)** | OCP with a free final time `tf`, solved by direct transcription and indirect shooting on a double integrator. |
| **[Free initial time](@ref tutorial-free-times-initial)** | OCP with a free initial time `t0`, solved by direct transcription and indirect shooting on a double integrator. |
| **[Free initial and final times](@ref tutorial-free-times-final-initial)** | OCP with both initial and final times as free optimization variables. |
| **[NLP manipulations](@ref tutorial-nlp)** | Low-level decomposition of `solve`: discretize the OCP, build an NLP model, solve it with an NLP solver, and rebuild the OCP solution. |
| **[Indirect simple shooting](@ref tutorial-indirect-simple-shooting)** | Implements the indirect simple shooting method on a simple example. |
| **[Goddard: direct, indirect](@ref tutorial-goddard)** | Solves the Goddard rocket problem (maximize altitude) with both direct and indirect methods, including singular and boundary arcs. |
| **[Linear–quadratic regulator](@ref tutorial-lqr)** | A simple LQR example illustrating how to solve a linear-quadratic OCP. |
| **[Minimal action](@ref tutorial-mam)** | Applies the Minimal Action Method (MAM) to find the most probable transition pathway between stable states in a stochastic system. |
| **[Model Predictive Control](@ref tutorial-mpc)** | Implements MPC for ship navigation in a sea current, with a real-time replanning loop. |
| **[Symbolic dynamics derivation](@ref tutorial-symbolics)** | Derives the equation of motion for the cart-pole symbolically and determines a periodic orbit. |

## Reproducibility

```@setup main
using Pkg
using InteractiveUtils
using Markdown

# Download links for the benchmark environment
function _downloads_toml(DIR)
    link_manifest = joinpath("assets", DIR, "Manifest.toml")
    link_project = joinpath("assets", DIR, "Project.toml")
    return Markdown.parse("""
    You can download the exact environment used to build this documentation:
    - 📦 [Project.toml]($link_project) - Package dependencies
    - 📋 [Manifest.toml]($link_manifest) - Complete dependency tree with versions
    """)
end
```

```@example main
_downloads_toml(".") # hide
```

```@raw html
<details style="margin-bottom: 0.5em; margin-top: 1em;"><summary>ℹ️ Version info</summary>
```

```@example main
versioninfo() # hide
```

```@raw html
</details>
```

```@raw html
<details style="margin-bottom: 0.5em;"><summary>📦 Package status</summary>
```

```@example main
Pkg.status() # hide
```

```@raw html
</details>
```

```@raw html
<details style="margin-bottom: 0.5em;"><summary>📚 Complete manifest</summary>
```

```@example main
Pkg.status(; mode = PKGMODE_MANIFEST) # hide
```

```@raw html
</details>
```
