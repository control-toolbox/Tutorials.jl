```@meta
Draft=false
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

const μ = 5.1658620912e12       # Constante gravitationnelle
const Tmax = 100.0              # Poussée max
const γmax = Tmax * 3600^2 / (2000 * 1e3)  # Contrôle max
const ε = 1e-2                  # Paramètre de régularisation

# État initial en coordonnées de Gauss 2D : (P, ex, ey, L)
x0 = [10000.0, 0.01, 0.0, 0.0]

# État final en coordonnées de Gauss 2D
xf = [42164.0, 0.0, 0.0, π]

T_min_100 = 13.4                # Temps minimal avec Tmax=100 (donné)
tf = 1.5 * T_min_100            # Temps final (à ajuster)
```

Dynamic functions in Gauss coordinates :

```@example orbit
asqrt(x; ε=1e-9) = sqrt(sqrt(x^2 + ε^2))  # sqrt lissée pour AD

function F0(x)
    P, ex, ey, L = x
    pdm = asqrt(P / μ)
    cl = cos(L)
    sl = sin(L)
    w = 1 + ex * cl + ey * sl
    F = zeros(eltype(x), 4)
    F[4] = w^2 / (P * pdm)     # dérivée de L (anomalie vraie)
    return F
end

function F1(x)
    # Direction de contrôle u₁ (dans le plan)
    P, ex, ey, L = x
    pdm = asqrt(P / μ)
    cl = cos(L)
    sl = sin(L)
    F = zeros(eltype(x), 4)
    F[2] = pdm * sl            # dérivée de ex
    F[3] = -pdm * cl           # dérivée de ey
    return F
end

function F2(x)
    # Direction de contrôle u₂ (dans le plan)
    P, ex, ey, L = x
    pdm = asqrt(P / μ)
    cl = cos(L)
    sl = sin(L)
    w = 1 + ex * cl + ey * sl
    F = zeros(eltype(x), 4)
    F[1] = pdm * 2 * P / w     # dérivée de P
    F[2] = pdm * (cl + (ex + cl) / w)
    F[3] = pdm * (sl + (ey + sl) / w)
    return F
end
```

Optimal control problem with regularization :

```@example orbit
function min_conso()
    @def ocp begin
        t ∈ [0, tf], time
        x = (P, ex, ey, L) ∈ R⁴, state
        u ∈ R², control

        # Conditions initiales et finales
        x(0) == x0
        x[1:3](tf) == xf[1:3]  # P, ex, ey à l'arrivée
        x[4](tf) == xf[4]      # anomalie vraie à l'arrivée

        # Dynamique
        ẋ(t) == F0(x(t)) + Tmax * (u₁(t) * F1(x(t)) + u₂(t) * F2(x(t)))

        norm = sqrt(u₁(t)^2 + u₂(t)^2)

        # Contrôle borné
        norm^2 ≤ γmax^2

        # Coût avec barrière logarithmique pour régularisation
        ∫(norm - ε * (log(norm) + log(γmax - norm))) → min
    end

    return ocp
end
nothing # hide
```

# Direct numerical solution

```@example orbit
using NLPModelsIpopt

u0 = [0.1, 0.0]
xguess(t) = x0 + (xf - x0) * t / tf

ocp = min_conso()
nlp_init = (state = xguess, control = u0)
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