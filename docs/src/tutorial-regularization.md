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


Tmax = 60                                  # Maximum thrust in Newtons
cTmax = 3600^2 / 1e6; T = Tmax * cTmax     # Conversion from Newtons to kg x Mm / h²
mass0 = 1500                               # Initial mass of the spacecraft
β = 1.42e-02                               # Engine specific impulsion
μ = 5165.8620912                           # Earth gravitation constant
P0 = 11.625                                # Initial semilatus rectum
ex0, ey0 = 0.75, 0                         # Initial eccentricity
hx0, hy0 = 6.12e-2, 0                      # Initial ascending node and inclination
L0 = π                                     # Initial longitude
Pf = 42.165                                # Final semilatus rectum
exf, eyf = 0, 0                            # Final eccentricity
hxf, hyf = 0, 0                            # Final ascending node and inclination

ε = 1e-1                             # Regularization parameter for logarithmic barrier
```

Dynamic functions in Gauss coordinates :

```@example orbit
asqrt(x; ε=1e-9) = sqrt(sqrt(x^2 + ε^2))  # sqrt lissée pour AD

function F0(x)
    P, ex, ey, hx, hy, L = x
    pdm = asqrt(P / μ)
    cl = cos(L)
    sl = sin(L)
    w = 1 + ex * cl + ey * sl
    F = zeros(eltype(x), 6) # Use eltype to allow overloading for AD
    F[6] = w^2 / (P * pdm)
    return F
end

function F1(x)
    P, ex, ey, hx, hy, L = x
    pdm = asqrt(P / μ)
    cl = cos(L)
    sl = sin(L)
    F = zeros(eltype(x), 6)
    F[2] = pdm *   sl
    F[3] = pdm * (-cl)
    return F
end

function F2(x)
    P, ex, ey, hx, hy, L = x
    pdm = asqrt(P / μ)
    cl = cos(L)
    sl = sin(L)
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
    cl = cos(L)
    sl = sin(L)
    w = 1 + ex * cl + ey * sl
    pdmw = pdm / w
    zz = hx * sl - hy * cl
    uh = (1 + hx^2 + hy^2) / 2
    F = zeros(eltype(x), 6)
    F[2] = pdmw * (-zz * ey)
    F[3] = pdmw *   zz * ex
    F[4] = pdmw *   uh * cl
    F[5] = pdmw *   uh * sl
    F[6] = pdmw *   zz
    return F
end
```

Initialisation with minimal time problem

```@example orbit
tf = 15                                      # Estimation of final time
Lf = 3π                                      # Estimation of final longitude
x0 = [P0, ex0, ey0, hx0, hy0, L0]            # Initial state
xf = [Pf, exf, eyf, hxf, hyf, Lf]            # Final state
x(t) = x0 + (xf - x0) * t / tf               # Linear interpolation
u = [0.1, 0.5, 0.]                        # Initial guess for the control
nlp_init = (state=x, control=u, variable=tf) # Initial guess for the NLP

function min_tf()
    @def ocp begin
        tf ∈ R, variable
        t ∈ [0, tf], time
        x = (P, ex, ey, hx, hy, L) ∈ R⁶, state
        u ∈ R³, control
        x(0) == x0
        x[1:5](tf) == xf[1:5]
        mass = mass0 - β * T * t
        ẋ(t) == F0(x(t)) + T / mass * (u₁(t) * F1(x(t)) + u₂(t) * F2(x(t)) + u₃(t) * F3(x(t)))
        u₁(t)^2 + u₂(t)^2 + u₃(t)^2 ≤ 1
        tf → min
    end
    return ocp
end
```

Optimal control problem with regularization :


```@example orbit
function min_conso()
    @def ocp begin
        t ∈ [0, tf], time
        x = (P, ex, ey, hx, hy, L) ∈ R⁶, state
        u ∈ R³, control
        x(0) == x0
        x[1:5](tf) == xf[1:5]
        mass = mass0 - β * T * t

        u_norm = sqrt(u₁(t)^2 + u₂(t)^2 + u₃(t)^2)

        # Dynamic
        ẋ(t) == F0(x(t)) + T / mass * (u₁(t) * F1(x(t)) + u₂(t) * F2(x(t)) + u₃(t) * F3(x(t)))
        1e-3 ≤ u_norm^2 ≤ 1
        u_g = max(1e-10, 1 - u_norm)


        # Regularization with logarithmic barrier
        ∫(u_norm - ε * (log(u_norm) + log(u_g))) → min
    end

    return ocp
end
```

# Direct numerical solution

```@example orbit
ocp1 = min_tf()  # Define the optimal control problem
sol = solve(ocp1; init=nlp_init, grid_size=100)

Tmax = 100.0              # Poussée max
T = Tmax * cTmax

tf = variable(sol)

ocp = min_conso()
nlp_init = (state = state(sol), control = control(sol))

nlp_sol = solve(ocp; init=nlp_init, grid_size=500)
plot(nlp_sol; control=:norm, size=(800, 300), layout=:group)
```

Plot in 3D :

```@example orbit
xsol = state(nlp_sol)
t = time_grid(nlp_sol)
N = size(t, 1)
nx = length(state(nlp_sol)(t[1]))  # nombre de variables d'état, par ex 6 ici

# Construire une matrice (nx x N) en évaluant xsol en chaque instant
X = zeros(nx, N)
for i in 1:N
    X[:, i] = xsol(t[i])
end

# Maintenant on peut accéder aux lignes comme prévu
P  = X[1, :]
ex = X[2, :]
ey = X[3, :]
hx = X[4, :]
hy = X[5, :]
L  = X[6, :]

cL = cos.(L)
sL = sin.(L)
w  = @. 1 + ex * cL + ey * sL
Z  = @. hx * sL - hy * cL
C  = @. 1 + hx^2 + hy^2

q1 = @. P *((1 + hx^2 - hy^2) * cL + 2 * hx * hy * sL) / (C * w)
q2 = @. P *((1 - hx^2 + hy^2) * sL + 2 * hx * hy * cL) / (C * w)
q3 = @. 2 * P * Z / (C * w)

plt1 = plot3d(1; xlim = (-60, 60), ylim = (-60, 60), zlim = (-5, 5), title = "Orbit transfer (direct)", legend=false)
@gif for i = 1:N
    push!(plt1, q1[i], q2[i], q3[i])
end every N ÷ min(N, 100)
```

# Indirect solution (Regularized shooting)
```@math
function ur(x, p)
    H1 = p' * F1(x)
    H2 = p' * F2(x)
    H3 = p' * F3(x)

    grad_H = [H1, H2, H3]
    normH = norm(grad_H)

    function φ(α)
        if α <= 1e-8 || α >= 1 - 1e-8
            return -Inf
        end
        return α * normH - ε * (log(α) + log(1 - α))
    end

    # Maximize φ(α) on [0,1]
    αs = range(1e-3, 1 - 1e-3, length=200)
    vals = φ.(αs)
    i = argmax(vals)
    αopt = αs[i]

    return αopt * grad_H / normH
end


fr = Flow(ocp, ur) # Regular flow (first version)

function shoot(ξ::Vector)
    tf = ξ[1]
    p0 = ξ[2:end]
    xf_, pf = fr(0, x0, p0, tf)

    s = [xf_[1:5] .- xf[1:5]; norm(p0)^2 - 1]

    return s
end

jshoot(ξ) = ForwardDiff.jacobian(shoot, ξ)
shoot!(s, ξ) = (s[:] = shoot(ξ); nothing)
jshoot!(js, ξ) = (js[:] = jshoot(ξ); nothing)

ξ0 = randn(6)
ξ0 ./= norm(ξ0)  # Normalize

# resolution of the shooting system with fsolve
bvp_sol = fsolve(shoot!, jshoot!, ξ0; show_trace=true)
println("Solution de la méthode de tir : ", bvp_sol)
```

# Conclusion

This tutorial shows how to solve a minimum-time orbital transfer with bounded thrust using a logarithmic barrier to regularize the problem. This approach facilitates the numerical solution by avoiding "bang-bang" control profiles.

The use of Gauss coordinates enables a simpler and more efficient modeling of orbital dynamics.

## References

Epenoy & Bertrand (CNES) "Optimal control and smoothing techniques..."

OptimalControl.jl Tutorial: Kepler minimum time

Cots Olivier, Contrôle optimal géométrique : méthodes homotopiques et applications (2012)