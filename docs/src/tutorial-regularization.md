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
using OrdinaryDiffEq
using Plots
using MINPACK
using ForwardDiff
using LinearAlgebra


γ_max = 0.1            # maximal thrust norm
const ε = 0.5           # Regularization parameter
const μ = 5.1658620912*1e12 # gravitational parameter
const P0 = 11625                                # Initial semilatus rectum
const ex0, ey0 = 0.75, 0                         # Initial eccentricity
const L0 = π                                     # Initial longitude
const Pf = 42165                                # Final semilatus rectum
const exf, eyf = 0, 0                            # Final eccentricity

tf = 30                                      # Estimation of final time
const Lf = 3π                                      # Estimation of final longitude
const x0 = [P0, ex0, ey0, L0]            # Initial state
const xf = [Pf, exf, eyf, Lf]            # Final state
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

        # Constrains
        x(0) == x0
        x[1:3](tf) == xf[1:3]  # P, ex, ey à l'arrivée
        x[4](tf) == xf[4]      # anomalie vraie à l'arrivée

        u_norm = sqrt(u₁(t)^2 + u₂(t)^2) + 1e-4

        # Control limit
        0.01 ≤ u_norm ≤ γmax - 0.01
        u_gamma = max(1e-10, γmax - u_norm)

        # Dynamic
        ẋ(t) == F0(x(t)) + Tmax * (u₁(t) * F1(x(t)) + u₂(t) * F2(x(t)))

        # Regularization with logarithmic barrier
        ∫(u_norm - ε * (log(u_norm) + log(u_gamma))) → min
    end

    return ocp
end
nothing # hide
```

Initialisation avec le problème à temps min

```@example orbit
function min_tf()
    @def ocp begin
        tf ∈ R, variable
        t ∈ [0, tf], time
        x = (P, ex, ey, L) ∈ R^4, state
        u ∈ R², control
        
        x(0) == x0
        x(tf) == xf
        0.05 ≤ tf
        
        ẋ(t) == F0(x(t)) + γ_max * (u₁(t) * F1(x(t)) + u₂(t) * F2(x(t)))
        u₁(t)^2 + u₂(t)^2 ≤ γ_max^2
        
        tf → min
    end

    return ocp
end
nothing # hide
```


# Direct numerical solution

```@example orbit
x(t) = x0 + (xf - x0) * t / tf               # Linear interpolation
u = [0.01, 0.05]                        # Initial guess for the control
init = (state=x, control=u, variable=tf) # Initial guess for the NLP
ocp_tf = min_tf()
sol = solve(ocp_tf; init=init, grid_size=100)

Tmax = 100.0              # Poussée max
γmax = Tmax * 3600^2 / 1e6  # Contrôle max

tf = variable(sol)

u0(t) = [0.05, 0.0]  # meilleur contrôle initial : norme non-nulle, mais pas au bord
function xguess(t)
    # Initialisation douce en évitant w ≈ 0
    α(t) = t / tf
    P = x0[1] + (xf[1] - x0[1]) * α(t)
    ex = (1 - α(t)) * x0[2]       # ex va vers 0 mais lentement
    ey = (1 - α(t)^2) * x0[3]     # idem
    L = α(t) * 2π + 0.1           # éviter cos(L) = -1
    return [P, ex, ey, L]
end

ocp = min_conso()
nlp_init = (state = xguess, control = u0)

nlp_sol = solve(ocp; init=nlp_init, grid_size=100, max_iter = 3000)
plot(nlp_sol)
```

# Indirect solution (Regularized shooting)

```@example orbit
function ur(x, p)
    H1 = p' * F1(x)
    H2 = p' * F2(x)

    normH = sqrt(H1^2 + H2^2)
    if normH < 1e-8
        return [0.0, 0.0]
    end
    γ = γmax
    u_norm = normH

    # optimal control with logarithm barrier
    u = [H1, H2] / u_norm
    return u
end


fr = Flow(ocp, ur) # Regular flow (first version)

function shoot(ξ::Vector)
    tf = ξ[1]
    p0 = ξ[2:end]
    xsol, psol = fr(0.0, x0, p0)
    xT = xsol[end]
    pT = psol[end]

    s = zeros(5)
    s[1:4] .= xT .- xf      # Conditions sur (P, ex, ey, L)
    s[5] = norm(p0)^2 - 1   # Normalisation du co-état
    return s
end

ξ = [T_min_100 * 1.5; randn(4)]
ξ[2:end] ./= norm(ξ[2:end])  # Normalisation initiale de p0
jshoot(ξ) = ForwardDiff.jacobian(shoot, ξ)
shoot!(s, ξ) = (s[:] = shoot(ξ); nothing)
jshoot!(js, ξ) = (js[:] = jshoot(ξ); nothing)
bvp_sol = fsolve(shoot!, jshoot!, ξ; show_trace=true); println(bvp_sol)
```

# Conclusion

This tutorial shows how to solve a minimum-time orbital transfer with bounded thrust using a logarithmic barrier to regularize the problem. This approach facilitates the numerical solution by avoiding "bang-bang" control profiles.

The use of Gauss coordinates enables a simpler and more efficient modeling of orbital dynamics.

## References

Epenoy & Bertrand (CNES) "Optimal control and smoothing techniques..."

OptimalControl.jl Tutorial: Kepler minimum time

Cots Olivier, Contrôle optimal géométrique : méthodes homotopiques et applications (2012)