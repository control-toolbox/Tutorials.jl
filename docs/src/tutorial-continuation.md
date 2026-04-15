# [Discrete continuation](@id tutorial-continuation)

```@meta
Draft = false
```

By using the warm start option, it is easy to implement a basic discrete continuation method, in which a sequence of problems is solved by using each solution as the initial guess for the next problem. This approach typically leads to faster and more reliable convergence than solving each problem with the same initial guess and is particularly useful for problems that require a good initial guess to converge.

## Continuation on parametric OCP

The most concise way to perform discrete continuation is to define a function that returns the optimal control problem for a given value of the continuation parameter, and then solve a sequence of such problems.
We illustrate this using a simple double integrator problem, where the fixed final time is gradually increased.

First we load the required packages:

```@example main-cont
using DataFrames
using OptimalControl
using NLPModelsIpopt
using Printf
using Plots
```

The `init` parameter of the `solve` function allows providing an initial guess. See the [initial guess documentation](@extref OptimalControl manual-initial-guess) for more details.

We write a function that returns the OCP for a given final time:

```@example main-cont
function problem(T)

    ocp = @def begin

        t ∈ [0, T], time
        x ∈ R², state
        u ∈ R, control

        q = x₁
        v = x₂

        q(0) == 0
        v(0) == 0
        q(T) == 1
        v(T) == 0
        ẋ(t) == [v(t), u(t)]

        ∫(u(t)^2) → min

    end

    return ocp
end
nothing # hide
```

Then we perform the continuation with a simple *for* loop, using each solution to initialize the next problem. We wrap the continuation in a function to avoid global variables.

```@example main-cont
function continuation_parametric(T_range; init=nothing, scheme=:midpoint, grid_size=200)
    data = DataFrame(T=Float64[], Objective=Float64[], Iterations=Int[])
    for T ∈ T_range
        ocp = problem(T)
        sol = solve(ocp; init=init, display=false, scheme=scheme, grid_size=grid_size)
        @assert successful(sol) "Solution failed for T=$T"
        init = sol
        push!(data, (T=T, Objective=objective(sol), Iterations=iterations(sol)))
    end
    return data
end

data = continuation_parametric(range(1, 2, length=5))
println(data)
```

We can visualize the evolution of the objective and the number of iterations with respect to the final time.

```@example main-cont
plt_obj = plot(data.T, data.Objective;
    seriestype=:scatter,
    title="Double integrator",
    label="Objective",
    xlabel="Final time T",
    ylabel="∫u(t)² dt")

plt_iter = plot(data.T, data.Iterations;
    seriestype=:scatter,
    label="Iterations",
    xlabel="Final time T",
    ylabel="Number of iterations")

layout = grid(2, 1, heights=[0.5, 0.5])
plot(plt_obj, plt_iter; layout=layout, size=(800, 600))
```

## Continuation on parameter

As a second example, we solve a Goddard problem with a decreasing maximum thrust. We define a function that returns the optimal control problem for a given value of `Tmax`, similar to the first example.

Let us first define the Goddard problem. Note that the formulation below illustrates all types of constraints, and the problem could be written more compactly.

```@example main-cont
# Parameters
r0 = 1
v0 = 0
m0 = 1
mf = 0.6
x0 = [r0, v0, m0]
vmax = 0.1

# Dynamics
function F0(x)
    # Uncontrolled dynamics: gravity and drag
    r, v, m = x
    D = Cd * v^2 * exp(-β*(r - 1))  # Aerodynamic drag
    return [ v, -D/m - 1/r^2, 0 ]   # [dr/dt, dv/dt, dm/dt]
end
function F1(x, Tmax)
    # Control dynamics: thrust contribution
    r, v, m = x
    return [ 0, Tmax/m, -b*Tmax ]   # [dr/dt, dv/dt, dm/dt] due to thrust
end

# Parameters for the dynamics
Cd = 310
β = 500
b = 2
Tmax_0 = 3
Tmax_f = 1

# Goddard problem function that takes Tmax as parameter
function goddard_problem(Tmax)
    ocp = @def begin

        tf ∈ R, variable
        t ∈ [0, tf], time
        x = (r, v, m) ∈ R^3, state
        u ∈ R, control

        0.01 ≤ tf ≤ Inf

        x(0) == x0
        m(tf) == mf
        r0 ≤ r(t) ≤ r0 + 0.1
        v0 ≤ v(t) ≤ vmax
        mf ≤ m(t) ≤ m0
        0 ≤ u(t) ≤ 1

        ẋ(t) == F0(x(t)) + u(t) * F1(x(t), Tmax)

        r(tf) → max

    end
    return ocp
end

# Solve the problem with a reference value of Tmax
sol0 = solve(goddard_problem(Tmax_0); display=false, scheme=:midpoint, grid_size=200)
@printf("Objective for reference solution: %.6f\n", objective(sol0))
```

Then, we perform the continuation on the maximal thrust. We wrap the continuation in a function that redefines the OCP at each step with the new `Tmax` value.

```@example main-cont
function continuation_goddard(Tmax_range; init=nothing, scheme=:midpoint, grid_size=200)
    data = DataFrame(Tmax=Float64[], Objective=Float64[], Iterations=Int[])
    sols = Vector{Any}()
    sol = init
    for Tmax_local ∈ Tmax_range
        ocp = goddard_problem(Tmax_local)
        sol = solve(ocp; init=sol, display=false, scheme=scheme, grid_size=grid_size)
        @assert successful(sol) "Solution failed for Tmax=$Tmax_local"
        push!(data, (Tmax=Tmax_local, Objective=objective(sol), Iterations=iterations(sol)))
        push!(sols, sol)
    end
    return data, sols
end

data, sols = continuation_goddard(range(Tmax_0, Tmax_f, length=9); init=sol0)
println(data)
```

We plot now the objective with respect to the maximal thrust, as well as the solutions for `Tmax=3`, `Tmax=2`, and `Tmax=1`. The time is normalized for the solution plots to compare trajectories with different final times.

```@example main-cont
using Plots.PlotMeasures # for leftmargin

plt_obj = plot(data.Tmax, data.Objective;
    seriestype=:scatter,
    title="Goddard problem",
    label="r(tf)",
    xlabel="Maximal thrust (Tmax)",
    ylabel="Maximal altitude r(tf)")

plt_iter = plot(data.Tmax, data.Iterations;
    seriestype=:scatter,
    label="Iterations",
    xlabel="Maximal thrust (Tmax)",
    ylabel="Number of iterations")

# Find indices closest to desired Tmax values
Tmax_values = [1.0, 2.0, 3.0]
indices = [argmin(abs.(data.Tmax .- T)) for T in Tmax_values]

layout_metrics = grid(2, 1, heights=[0.5, 0.5])
plot(plt_obj, plt_iter; layout=layout_metrics, size=(800, 600), leftmargin=5mm)
```

We now plot the solutions for the selected `Tmax` values, using normalized time to compare trajectories with different final times.

```@example main-cont
plt_sol = plot()
for (i, idx) in enumerate(indices)
    plot!(plt_sol, sols[idx]; label="Tmax=$(round(data.Tmax[idx], digits=2))", time=:normalize, color=i)
end
plot(plt_sol; size=(800, 800), leftmargin=5mm)
```

### Animation

The animation below shows three rockets launching simultaneously with different maximum thrust values (Tmax = 3, 2, 1). The rockets start at the same time (t=0) but reach their final altitudes at different times according to their optimal trajectories. Each rocket displays:

- **Altitude** (above the rocket)
- **Velocity gauge** (horizontal bar, green→red)
- **Fuel gauge** (vertical bar, showing remaining fuel m(t) - mf)
- **Flame** (visible when thrust u(t) > 0)

```@setup main-cont
# Extract data for the 3 selected solutions (Tmax ≈ 3, 2, 1)
selected_sols = [sols[idx] for idx in indices]
selected_Tmax = [round(data.Tmax[idx], digits=2) for idx in indices]

# Extract time grids and state/control data for each solution
solutions_data = []
for (sol, Tmax_val) in zip(selected_sols, selected_Tmax)
    t_grid = time_grid(sol)
    X = state(sol).(t_grid)
    u = control(sol).(t_grid)
    
    # Extract r, v, m from state
    r = [x[1] for x in X]
    v = [x[2] for x in X]
    m = [x[3] for x in X]
    
    push!(solutions_data, (t=t_grid, r=r, v=v, m=m, u=u))
end

# Serialize to JSON for JavaScript injection
json_t1 = "[" * join(string.(solutions_data[1].t), ",") * "]"
json_r1 = "[" * join(string.(solutions_data[1].r), ",") * "]"
json_v1 = "[" * join(string.(solutions_data[1].v), ",") * "]"
json_m1 = "[" * join(string.(solutions_data[1].m), ",") * "]"
json_u1 = "[" * join(string.(solutions_data[1].u), ",") * "]"

json_t2 = "[" * join(string.(solutions_data[2].t), ",") * "]"
json_r2 = "[" * join(string.(solutions_data[2].r), ",") * "]"
json_v2 = "[" * join(string.(solutions_data[2].v), ",") * "]"
json_m2 = "[" * join(string.(solutions_data[2].m), ",") * "]"
json_u2 = "[" * join(string.(solutions_data[2].u), ",") * "]"

json_t3 = "[" * join(string.(solutions_data[3].t), ",") * "]"
json_r3 = "[" * join(string.(solutions_data[3].r), ",") * "]"
json_v3 = "[" * join(string.(solutions_data[3].v), ",") * "]"
json_m3 = "[" * join(string.(solutions_data[3].m), ",") * "]"
json_u3 = "[" * join(string.(solutions_data[3].u), ",") * "]"

# Inject real Tmax values for labels
json_Tmax1 = string(selected_Tmax[1])
json_Tmax2 = string(selected_Tmax[2])
json_Tmax3 = string(selected_Tmax[3])

# Inject m0 and mf values
json_m0 = string(m0)
json_mf = string(mf)

# Inject r0 value
json_r0 = string(r0)

# Calculate global max time for animation loop
t_max_global = max(solutions_data[1].t[end], solutions_data[2].t[end], solutions_data[3].t[end])

# Define RawHTML wrapper for Documenter
struct RawHTML
    raw::String
end
Base.show(io::IO, ::MIME"text/html", h::RawHTML) = print(io, h.raw)

html_anim = """
<div style="display: flex; justify-content: center; margin: 20px 0;">
    <canvas id="goddardCanvas" width="900" height="650"
            style="border:1px solid #ddd; border-radius: 8px; max-width: 100%;">
    </canvas>
</div>

<script>
(function() {
    // Data injected from Julia
    const t1 = $json_t1, r1 = $json_r1, v1 = $json_v1, m1 = $json_m1, u1 = $json_u1;
    const t2 = $json_t2, r2 = $json_r2, v2 = $json_v2, m2 = $json_m2, u2 = $json_u2;
    const t3 = $json_t3, r3 = $json_r3, v3 = $json_v3, m3 = $json_m3, u3 = $json_u3;
    
    const Tmax1 = $json_Tmax1, Tmax2 = $json_Tmax2, Tmax3 = $json_Tmax3;
    
    const t_max_global = $t_max_global;
    const m0 = $json_m0, mf = $json_mf;
    const r0 = $json_r0;
    
    const canvas = document.getElementById('goddardCanvas');
    const ctx = canvas.getContext('2d');
    
    // Detect dark/light mode
    const isDark = window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches;
    
    // Color palettes
    const colors = isDark ? {
        canvas_bg: '#1e1e2e',
        ground: '#4a4a4a',
        gauge_bg: '#3a3a3a',
        gauge_high: '#2ecc71',
        gauge_mid: '#f39c12',
        gauge_low: '#e74c3c',
        text: '#ecf0f1',
        text_secondary: '#bdc3c7',
        flame_outer: '#e67e22',
        flame_inner: '#f1c40f',
        window_fill: 'rgba(255,255,255,0.3)',
        window_stroke: 'rgba(255,255,255,0.6)',
        progress: '#4063D8'
    } : {
        canvas_bg: '#fafafa',
        ground: '#ccc',
        gauge_bg: '#ecf0f1',
        gauge_high: '#2ecc71',
        gauge_mid: '#f39c12',
        gauge_low: '#e74c3c',
        text: '#2c3e50',
        text_secondary: '#7f8c8d',
        flame_outer: '#e67e22',
        flame_inner: '#f1c40f',
        window_fill: 'rgba(255,255,255,0.55)',
        window_stroke: 'rgba(255,255,255,0.85)',
        progress: '#4063D8'
    };
    
    // Animation duration in seconds (independent of real problem time)
    const animation_duration = 10.0;
    
    const zone_width = canvas.width / 3;
    const top_margin = 50;   // Space for time label at top
    const bottom_margin = 120; // Space for ground line + v gauge + Tmax label
    const baseline = canvas.height - bottom_margin;
    
    // Scale for altitude (r goes from r0 to ~r0+0.1)
    const r_min = r0;
    const r_max_all = Math.max(...r1, ...r2, ...r3);
    const r_range = r_max_all - r_min;
    const altitude_scale = (baseline - top_margin) / r_range;
    
    // Scale for velocity (v goes from 0 to ~0.1)
    const v_max_all = Math.max(...v1, ...v2, ...v3);
    
    let start_time = null;

    function interpolate(t, data_t, data_y) {
        // Find interval
        let i = 0;
        while(i < data_t.length - 1 && data_t[i + 1] < t) { i++; }
        
        if (i >= data_t.length - 1) {
            return data_y[data_y.length - 1];  // Return final value if past end
        }
        
        const dt = data_t[i+1] - data_t[i];
        const alpha = (dt > 0) ? (t - data_t[i]) / dt : 0;
        return data_y[i] + alpha * (data_y[i+1] - data_y[i]);
    }

    function drawRocket(zone_idx, t, r, v, m, u, color, Tmax) {
        const zone_x = zone_idx * zone_width;
        const center_x = zone_x + zone_width / 2;
        
        // Rocket position
        const rocket_y = baseline - altitude_scale * (r - r_min);
        
        // Draw ground line
        ctx.beginPath();
        ctx.moveTo(zone_x + 20, baseline + 20);
        ctx.lineTo(zone_x + zone_width - 20, baseline + 20);
        ctx.lineWidth = 2;
        ctx.strokeStyle = colors.ground;
        ctx.stroke();
        
        // Draw fuel gauge (vertical bar to the right of rocket)
        const fuel_x = center_x + 80;
        const fuel_height = 200;
        const fuel_y = baseline - fuel_height;
        const fuel_level = (m - mf) / (m0 - mf);  // 0 to 1
        
        // Fuel gauge background
        ctx.fillStyle = colors.gauge_bg;
        ctx.fillRect(fuel_x - 10, fuel_y, 20, fuel_height);
        
        // Fuel gauge fill
        const fuel_fill_height = fuel_level * fuel_height;
        ctx.fillStyle = fuel_level > 0.5 ? colors.gauge_high : (fuel_level > 0.2 ? colors.gauge_mid : colors.gauge_low);
        ctx.fillRect(fuel_x - 10, fuel_y + (fuel_height - fuel_fill_height), 20, fuel_fill_height);
        
        // Fuel gauge label
        ctx.fillStyle = colors.text_secondary;
        ctx.font = '14px Arial';
        ctx.textAlign = 'center';
        ctx.fillText('Fuel', fuel_x, fuel_y - 8);
        
        // Draw rocket: nose cone + body + fins + window
        const rocket_bottom = rocket_y + 18;
        ctx.fillStyle = color;

        // Nose cone
        ctx.beginPath();
        ctx.moveTo(center_x, rocket_y - 38);
        ctx.lineTo(center_x - 12, rocket_y - 10);
        ctx.lineTo(center_x + 12, rocket_y - 10);
        ctx.closePath();
        ctx.fill();

        // Body
        ctx.fillRect(center_x - 12, rocket_y - 10, 24, 28);

        // Left fin
        ctx.beginPath();
        ctx.moveTo(center_x - 12, rocket_y + 6);
        ctx.lineTo(center_x - 24, rocket_y + 22);
        ctx.lineTo(center_x - 12, rocket_y + 18);
        ctx.closePath();
        ctx.fill();

        // Right fin
        ctx.beginPath();
        ctx.moveTo(center_x + 12, rocket_y + 6);
        ctx.lineTo(center_x + 24, rocket_y + 22);
        ctx.lineTo(center_x + 12, rocket_y + 18);
        ctx.closePath();
        ctx.fill();

        // Window
        ctx.beginPath();
        ctx.arc(center_x, rocket_y + 2, 5, 0, 2 * Math.PI);
        ctx.fillStyle = colors.window_fill;
        ctx.fill();
        ctx.strokeStyle = colors.window_stroke;
        ctx.lineWidth = 1.5;
        ctx.stroke();

        // Draw flame if thrust > 0.1
        if (u > 0.1) {
            const flame_size = 10 + 20 * u;
            ctx.fillStyle = colors.flame_outer;
            ctx.beginPath();
            ctx.moveTo(center_x - 10, rocket_bottom);
            ctx.lineTo(center_x, rocket_bottom + flame_size);
            ctx.lineTo(center_x + 10, rocket_bottom);
            ctx.closePath();
            ctx.fill();

            // Inner flame
            ctx.fillStyle = colors.flame_inner;
            ctx.beginPath();
            ctx.moveTo(center_x - 5, rocket_bottom);
            ctx.lineTo(center_x, rocket_bottom + flame_size * 0.6);
            ctx.lineTo(center_x + 5, rocket_bottom);
            ctx.closePath();
            ctx.fill();
        }
        
        // Draw altitude text
        ctx.fillStyle = colors.text;
        ctx.font = 'bold 14px Arial';
        ctx.textAlign = 'center';
        ctx.fillText('r = ' + r.toFixed(4), center_x, rocket_y - 40);
        
        // Draw velocity gauge (horizontal bar)
        const v_gauge_width = 60;
        const v_gauge_height = 8;
        const v_gauge_x = center_x - v_gauge_width / 2;
        const v_gauge_y = rocket_y + 50;
        const v_level = Math.abs(v) / v_max_all;
        
        // Velocity gauge background
        ctx.fillStyle = colors.gauge_bg;
        ctx.fillRect(v_gauge_x, v_gauge_y, v_gauge_width, v_gauge_height);
        
        // Velocity gauge fill
        ctx.fillStyle = v_level > 0.7 ? colors.gauge_low : (v_level > 0.4 ? colors.gauge_mid : colors.gauge_high);
        ctx.fillRect(v_gauge_x, v_gauge_y, v_level * v_gauge_width, v_gauge_height);
        
        // Velocity label
        ctx.fillStyle = colors.text_secondary;
        ctx.font = '12px Arial';
        ctx.fillText('v', center_x, v_gauge_y + 20);
        
        // Draw Tmax label at bottom
        ctx.fillStyle = colors.text;
        ctx.font = 'bold 16px Arial';
        ctx.fillText('Tmax = ' + Tmax, center_x, baseline + 90);
    }

    function draw(time) {
        if (!start_time) start_time = time;
        
        const elapsed = (time - start_time) / 1000.0;
        const progress = (elapsed % animation_duration) / animation_duration;
        const t_anim = progress * t_max_global;
        
        ctx.fillStyle = colors.canvas_bg;
        ctx.fillRect(0, 0, canvas.width, canvas.height);
        
        // Draw time indicator at top
        ctx.fillStyle = colors.text;
        ctx.font = 'bold 18px Arial';
        ctx.textAlign = 'center';
        ctx.fillText('t = ' + t_anim.toFixed(3) + ' s', canvas.width / 2, 30);
        
        // Interpolate and draw each rocket
        const r1_cur = interpolate(t_anim, t1, r1);
        const v1_cur = interpolate(t_anim, t1, v1);
        const m1_cur = interpolate(t_anim, t1, m1);
        const u1_cur = interpolate(t_anim, t1, u1);
        drawRocket(0, t_anim, r1_cur, v1_cur, m1_cur, u1_cur, '#CB3C33', Tmax1);
        
        const r2_cur = interpolate(t_anim, t2, r2);
        const v2_cur = interpolate(t_anim, t2, v2);
        const m2_cur = interpolate(t_anim, t2, m2);
        const u2_cur = interpolate(t_anim, t2, u2);
        drawRocket(1, t_anim, r2_cur, v2_cur, m2_cur, u2_cur, '#389826', Tmax2);
        
        const r3_cur = interpolate(t_anim, t3, r3);
        const v3_cur = interpolate(t_anim, t3, v3);
        const m3_cur = interpolate(t_anim, t3, m3);
        const u3_cur = interpolate(t_anim, t3, u3);
        drawRocket(2, t_anim, r3_cur, v3_cur, m3_cur, u3_cur, '#9558B2', Tmax3);
        
        // Draw progress bar at bottom
        const bar_height = 4;
        ctx.fillStyle = colors.gauge_bg;
        ctx.fillRect(0, canvas.height - bar_height, canvas.width, bar_height);
        ctx.fillStyle = colors.progress;
        ctx.fillRect(0, canvas.height - bar_height, canvas.width * progress, bar_height);
        
        requestAnimationFrame(draw);
    }
    requestAnimationFrame(draw);
})();
</script>
"""
```

```@example main-cont
RawHTML(html_anim)  # hide
```

## Best practices and limitations

The examples above illustrate a basic discrete continuation method. Here are some important considerations:

- **Step size**: The continuation step size should be chosen carefully. Too large a step may lead to convergence issues, while too small a step may be inefficient. Adaptive step size strategies can improve robustness but are not covered in this tutorial.

- **Convergence monitoring**: Always check that each solution converges successfully using `successful(sol)`. If a solution fails, consider reducing the step size or providing a better initial guess.

- **Advanced methods**: For more challenging problems, homotopic continuation methods with path following algorithms can be used. These methods adapt the continuation parameter automatically and can handle bifurcations. Such advanced techniques are beyond the scope of this tutorial.

- **Alternative approaches**: In some cases, it may be more efficient to directly parameterize the continuation parameter within the optimal control problem rather than solving a sequence of problems. OptimalControl plans to add this feature in a future release.
