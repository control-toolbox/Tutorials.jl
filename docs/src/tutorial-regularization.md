```@meta
    draft=false
```
## Regularization method application :

Regularization in optimal control helps to smooth discontinuous or “bang-bang” control laws, making the problem more numerically tractable. It improves the convergence of solvers by ensuring differentiability. Regularization also avoids issues like singularities or ill-conditioned Jacobians in shooting methods. A common technique is to penalize controls near their bounds using barrier or penalty functions. This allows gradually approaching the true optimal solution while maintaining stability during computation.

## Regularization of an Orbital Transfer with Logarithmic Barrier

# Introduction

This tutorial demonstrates how to regularize a bounded-thrust orbital transfer optimal control problem, typically of "bang-bang" type, using a logarithmic barrier on the control. We will use Gauss coordinates, which simplify the equations of orbital motion.

The problem will be solved using two approaches:

    - a direct method (NLP formulation + Ipopt solver),

    - an indirect method (shooting with a regularized Hamiltonian).

Problem data and constants :

```@example orbit
using OptimalControl
using NLPModelsIpopt
using BenchmarkTools
using DataFrames
using Plots
using Printf
using LinearAlgebra

const Tmax = 60                                # Maximum thrust (N)
const cTmax = 3600^2 / 1e6; const T = Tmax * cTmax  # Conversion N → kg*Mm/h²
const mass0 = 1500                             # Initial mass (kg)
const β = 1.42e-2                              # Specific impulse
const μ = 5165.8620912                         # Gravitational constant

const P0 = 11.625
const ex0, ey0 = 0.75, 0
const hx0, hy0 = 6.12e-2, 0
const L0 = π
const Pf = 42.165
const exf, eyf = 0.0, 0.0
const hxf, hyf = 0.0, 0.0
const ε = 1e-2                                 # Regularization parameter

x0 = [P0, ex0, ey0, hx0, hy0, L0]
xf = [Pf, exf, eyf, hxf, hyf]
```

Dynamic functions in Gauss coordinates :

```@example orbit
asqrt(x; ε=1e-9) = sqrt(sqrt(x^2 + ε^2))

function F0(x)
    P, ex, ey, hx, hy, L = x
    pdm = asqrt(P / μ)
    cl = cos(L); sl = sin(L)
    w = 1 + ex * cl + ey * sl
    F = zeros(eltype(x), 6)
    F[6] = w^2 / (P * pdm)
    return F
end

function F1(x)
    P, ex, ey, hx, hy, L = x
    pdm = asqrt(P / μ)
    cl = cos(L); sl = sin(L)
    F = zeros(eltype(x), 6)
    F[2] = pdm * sl
    F[3] = pdm * (-cl)
    return F
end

function F2(x)
    P, ex, ey, hx, hy, L = x
    pdm = asqrt(P / μ)
    cl = cos(L); sl = sin(L)
    w = 1 + ex * cl + ey * sl
    F = zeros(eltype(x), 6)
    F[1] = pdm * 2 * P / w
    F[2] = pdm * (cl + (ex + cl) / w)
    F[3] = pdm * (sl + (ey + sl) / w)
    return F
end

function F3(x)
    P, ex, ey, hx, hy, L = x
    pdm = asqrt(P / μ)
    cl = cos(L); sl = sin(L)
    w = 1 + ex * cl + ey * sl
    pdmw = pdm / w
    zz = hx * sl - hy * cl
    uh = (1 + hx^2 + hy^2) / 2
    F = zeros(eltype(x), 6)
    F[2] = pdmw * (-zz * ey)
    F[3] = pdmw * (zz * ex)
    F[4] = pdmw * uh * cl
    F[5] = pdmw * uh * sl
    F[6] = pdmw * zz
    return F
end
```

Optimal control problem with regularization :

```@example orbit
@def ocp begin
    tf ∈ R, variable
    t ∈ [0, tf], time
    x = (P, ex, ey, hx, hy, L) ∈ R⁶, state
    u ∈ R³, control
    x(0) == x0
    x[1:5](tf) == xf[1:5]
    mass = mass0 - β * T * t
    ẋ(t) == F0(x(t)) + (T / mass) * (u₁(t) * F1(x(t)) + u₂(t) * F2(x(t)) + u₃(t) * F3(x(t)))
    u₁(t)^2 + u₂(t)^2 + u₃(t)^2 < 1
    tf → min + ε * ∫(log(1 - (u₁(t)^2 + u₂(t)^2 + u₃(t)^2)))
end
```

# Direct numerical solution

```@example orbit
using NLPModelsIpopt

u0 = [0.1, 0.5, 0.0]
xguess(t) = x0 + (xf - x0) * t / 15
nlp_init = (state=xguess, control=u0, variable=15.0)
nlp_sol = solve(ocp; init=nlp_init, grid_size=100)
plot(nlp_sol)
```

# Indirect solution (Regularized shooting)



# Conclusion

This tutorial shows how to solve a minimum-time orbital transfer with bounded thrust using a logarithmic barrier to regularize the problem. This approach facilitates the numerical solution by avoiding "bang-bang" control profiles.

The use of Gauss coordinates enables a simpler and more efficient modeling of orbital dynamics.

## References

Epenoy & Bertrand (CNES) "Optimal control and smoothing techniques..."

OptimalControl.jl Tutorial: Kepler minimum time

Cots Olivier, Contrôle optimal géométrique : méthodes homotopiques et applications (2012)