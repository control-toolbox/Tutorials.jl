# [Optimal control problem with free final time](@id tutorial-free-times-final)

```@meta
Draft = true
```

In this tutorial, we explore an optimal control problem with free final time `tf`.

Here is the required packages for the tutorial:

```@example main-free-final
using LinearAlgebra: norm
using NLPModelsIpopt
using NonlinearSolve
using OptimalControl
using OrdinaryDiffEq  # to get the Flow function from OptimalControl
using Plots
using Printf
```

## Definition of the problem

```@example main-free-final
t0 = 0              # initial time
x0 = [0, 0]         # initial point
xf_target = [1, 0]  # final point

@def ocp begin
    tf ∈ R, variable
    t ∈ [t0, tf], time
    x ∈ R², state
    u ∈ R, control
    -1 ≤ u(t) ≤ 1
    x(t0) == x0
    x(tf) == xf_target
    0.05 ≤ tf ≤ Inf
    ẋ(t) == [x₂(t), u(t)]
    tf → min
end
nothing # hide
```

In the definition, the line that shows that the problem has a **free** final time is the following:

```julia
tf ∈ R, variable
```

The variable `tf` is indeed the final time according to:

```julia
t ∈ [t0, tf], time
```

## Direct resolution

You may need to add this type of constrain in your problem definition :

To help the direct solver to converge, we have added a constraint on the variable `tf`:

```julia
0.05 ≤ tf ≤ Inf
```

We now solve the problem using a direct method, with automatic treatment of the free final time.

```@example main-free-final
sol = solve(ocp; grid_size=100)
```

And plot the solution.

```@example main-free-final
plt = plot(sol; label="direct", size=(800, 800))
```

## Verification of results

!!! note "Mathematical computations"

    ```@raw html
    <details>
    <summary> Click to show/hide mathematical computations.</summary>
    ```

    Here is the theoretical part. The pseudo-Hamiltonian is:

    ```math
        H(x, p, u) = p_1 x_2 + p_2 u.
    ```

    Pontryagin's theorem gives:

    ```math
        \begin{aligned}
            \dot{p}_1(t) &= 0, \\
            \dot{p}_2(t) &= -p_1(t). \\
        \end{aligned}
    ```

    Hence, $p_1(t) = \mathrm{cst} = \alpha$ and $p_2(t) = -\alpha t + \beta$, $\beta = p_2(0)$. According to Pontryagin's theorem, we have:

     ```math
        \begin{aligned}
            p_1(t) &= 1, \\
            p_2(t) &= 1-t. \\
        \end{aligned}
    ```   

    Besides, the pseudo-Hamiltonian satisfies:
    ```math
        H(x(t_f), p(t_f), u(t_f)) = -p° = 1.
    ```

    For this problem, we can prove that the optimal control $u$ satisfies:

    ```math
        u(t) = \left\{
        \begin{aligned}
             1 & \quad\text{if}\quad t \in [0, t_1], \\
            -1 & \quad\text{if}\quad t \in (t_1, t_f], \\
        \end{aligned}
        \right.
    ```
    
    where $t_1$ is a constant to determined, as $t_f$. To do so, we have to compute the state $x$:

    On t ∈ [0,t1] :
    ```math
        \begin{aligned}
            \dot{x}_1(t) &= x_2(t) \\
            \dot{x}_2(t) &= 1 \\
        \end{aligned}
    ```

    ```math
        \begin{aligned}
            x_2(t) &= t \\
            x_1(t) &= (1/2)*t^2
        \end{aligned}

    ```

    When t = t1 :
    ```math
        \begin{aligned}
            x_2(t1) &= t_1 \\
            x_1(t1) &= (1/2)*t_1^2
        \end{aligned}
    ```

    On $t \in [t_1,t_f]$
    ```math
        \begin{aligned}
            x_2(t) &= -t + 2t_1 \\
            x_1(t) &= -\frac{1}{2}t^2 + 2t_1 t + C_2 \\
        \end{aligned}
    ```
    ```math
        \begin{aligned}
            x_1(t_1) &= -\frac{1}{2}t_1^2 + 2t_1^2 + C_2 &= \frac{1}{2}t_1^2 \\ 
            C_2 &= -t_1^2 \\
        \end{aligned}
    ```

    ```math
        \begin{aligned}
            x_1(t) &= -\frac{1}{2}t^2 + 2t_1 t - t_1^2
        \end{aligned}
    ```
    Finally you can solve the terminal conditions :
    ```math
        \begin{aligned}
            x_1(t_f) &= 1, \quad x_2(t_f) &= 0 \\
        \end{aligned}
    ```
    ```math
        \begin{aligned}
            x_2(t_f) &= -t_f + 2t_1 &= 0 \\
            t_f &= 2t_1 \\
        \end{aligned}
    ```
    ```math
        \begin{aligned}
            x_1(t_f)& = -\frac{1}{2}t_f^2 + 2t_1 t_f - t_1^2 &= 1 \\
        \end{aligned}
    ```
    and 
    ```math
        \begin{aligned}
            -2t_1^2 + 4t_1^2 - t_1^2 &= t_1^2 &= 1 \\
        \end{aligned}
    ```
    gets us
    ```math
        \begin{aligned}
            t_1 &= 1, \quad t_f &= 2
        \end{aligned}
    ```
    To sum up we find the following solutions :
    ```math

    x_1(t) = \begin{cases}
    \frac{1}{2}t^2 & \text{si } t \in [0, 1) \\
    -\frac{1}{2}t^2 + 2t - 1 & \text{si } t \in [1, 2]
    \end{cases}
    \qquad
    x_2(t) = \begin{cases}
    t & \text{si } t \in [0, 1) \\
    2 - t & \text{si } t \in [1, 2]
    \end{cases}
    ```
    ```math
    p_1(t) = 1, \quad p_2(t) = 1-t, \quad p^0=-1
    ```
    ```math
    u(t) = \begin{cases}
    1 & \text{si } t \in [0, 1) \\
    -1 & \text{si } t \in [1, 2]
    \end{cases}
    ```
    ```@raw html
    </details>
    ```

Now we can compare the results found with the direct method with the theoretical analysis:

```@example main-free-final
tf = variable(sol)
u = control(sol)
p = costate(sol)
x = state(sol)
p° = -1
H(t) = p(t)[1]*x(t)[2] + p(t)[2]*u(t) 

@printf("H(tf) = %.3f\n", H(tf))
@printf("x(tf) = [%.3f, %.3f]\n", x(tf)[1], x(tf)[2])
@printf("p(tf) = [%.3f, %.3f]\n", p(tf)[1], p(tf)[2])
```

The numerical results closely match the theoretical predictions: the final state $x(t_f)=[1,0]$ is exactly satisfied.
The costate and Hamiltonian values at final time show a small deviation (≈ 0.01), likely due to numerical precision.
Overall, the direct method confirms the theoretical analysis with excellent accuracy.

We can analyze the influence of using different discretization sizes (`grid_size`), and observed the following results for the optimal $t_f$:

```@example main-free-final
for N in [20, 50, 100, 200]
    solN = solve(ocp; grid_size=N, display=false)
    @printf("grid_size = %3d → tf = %.5f\n", N, objective(solN))
end
```

This example shows that problems with a free final time can be sensitive to discretization. A small grid may lead to suboptimal or slightly inaccurate results.

## Indirect method

We first define the pseudo-Hamiltonian and the switching function needed to define the shooting function.

```@example main-free-final
# Pseudo-Hamitlonian
H(x, p, u) = p[1]*x[2] + p[2]*u

# Hamiltonian lift
F1(x) = [0, 1]
H1 = Lift(F1)
nothing # hide
```

Then, we define the different flows, associated to control laws $u(t)=+1$ and $u(t)=-1$.

```@example main-free-final
const u_pos = 1
const u_neg = -1

f_pos = Flow(ocp, (x, p, tf) -> u_pos) # it depends on tf since tf is the variable
f_neg = Flow(ocp, (x, p, tf) -> u_neg)
nothing # hide 
```

We can now define the shooting function following the structure given by the direct method.

```@example main-free-final
function shoot!(s, p0, t1, tf)

    # the flows
    x1, p1 = f_pos(t0, x0, p0, t1)
    xf, pf = f_neg(t1, x1, p1, tf)

    # final control
    uf = -1

    # shooting conditions
    s[1:2] = xf .- xf_target    # reach the target
    s[3]   = H(xf, pf, uf) - 1  # final condition on the pseudo-Hamiltonian
    s[4]   = H1(x1, p1)         # switching condition

end
nothing # hide
```

Before solving our problem we must find a good initial guess to help the convergence of the algorithm.

```@example main-free-final
t = time_grid(sol)
x = state(sol)
p = costate(sol)
φ(t) = H1(x(t), p(t))  # switching function

p0 = p(t0)
t1 = t[argmin(abs.(φ.(t)))]
tf = t[end]

println("p0 = ", p0)
println("t1 = ", t1)
println("tf = ", tf)

# Norm of the shooting function at initial guess
s = zeros(4)
shoot!(s, p0, t1, tf)
println("\nNorm of the shooting function: ‖s‖ = ", norm(s), "\n")
```

And we finally can solve the problem using an indirect method.

```@example main-free-final
# Aggregated function
nle!  = (s, ξ, λ) -> shoot!(s, ξ[1:2], ξ[3], ξ[4])

# Nonlinear problem
ξ_guess = [p0..., t1, tf] # Initial guess
prob = NonlinearProblem(nle!, ξ_guess)

# Resolution
indirect_sol = solve(prob; show_trace=Val(true), abstol=1e-8, reltol=1e-8)
```

We can now plot and compare with the direct method.

```@example main-free-final
# Data
p0 = indirect_sol.u[1:2]
t1 = indirect_sol.u[3]
tf = indirect_sol.u[4]

# Compute the optimal solution
f = f_pos * (t1, f_neg) # concatenate the flows
flow_sol = f((t0, tf), x0, p0; saveat=range(t0, tf, 100))

# Plot the solution
plot!(plt, flow_sol; label="indirect", color=:red)
```
