# [constraints-at-intermediate-times](@id tutorial-for-constraints-at-intermediate-times)
In this tutorial, we will explore an example for constraints at intermediate times
We want to showcase use cases where we have some intermediate times with constraints associated to each one:


## General Problem Formulation

The optimal control problem is defined as:
```math
\min_{q,u} \int_{t_0}^{t_f} f_0(x(t), u(t)) dt 
```

Subject to the dynamic constraint:
```math
\dot{x}(t) = f(x(t), u(t)), \quad u(t) \in U
```

With boundary conditions:
```math
c_s(x(t_s)) = 0, \quad s = 0, ..., N
```

## Step-by-Step Derivation

### Step 1: Time Discretization
With discretization points:
```math
\text{avec} \quad t_0 < ... < t_k < ... < t_f = t_N
```

For N = 2, we define:
```math
 q_1(s) = x(t) \quad \text{pour } t \in [t_0, t_1]
 q_2(s) = x(t_0^s) \quad \text{pour } s \in [t_1, t_2]

```

### Step 2: Discretized System


Therefore:
```math
\Rightarrow  q_k(s) = x(\varphi_k(s)) \quad \text{for } s \in [0,1] \text{ and with } \varphi_k(s) = t_0 + (t_k - t_(k-1)s
```

Where:
```math
\varphi_1(s) = \varphi(t_0, t_1, s)
```
```math
\varphi_2(s) = \varphi(t_1, t_2, s)
```

### Interpolation Functions
The interpolation functions are defined as:
```math
\Rightarrow \psi(a,b,s) = a + (b-a)s
```


## Final Problem Definition

**minimize:**
```math
\min_{q,v} \sum_{k=1}^{N} \Delta t_k \int_0^1 f_0(q_k(s), v_k(s)) ds 
```

**subject to:**
```math
q_k(s) = \Delta t_k  f(q_k(s), u_k(s)) ds \quad \forall s \in [0,1]
```

**with constraints:**
```math
c_s(q_1(0)) = 0
```
```math
c_s(q_k(1)) = 0 \quad \forall k \in [0,N]
```
```math
q_k(1) = q_{k+1}(0) \quad \forall k \in [0,N-1]
```
## Example to modelize a trein passing through multiple stations 

Necessary imports:

```@example main-cit
using OptimalControl
using NLPModelsIpopt
using Plots
using Plots.PlotMeasures # for leftmargin, bottommargin

# Tram system parameters
max_acceleration = 2.0  # m/s^2
max_deceleration = -3.0 # m/s^2
max_velocity = 20.0     # m/s (72 km/h)

# Station positions along the track (in meters)
station_0 = 0.0      # Starting station
station_1 = 800.0    # First intermediate station
station_2 = 1600.0   # Second intermediate station  
station_3 = 2400.0   # Final destination

# Time windows for each arc (in seconds)
t0 = 0.0
t1 = 60.0   # Time to reach first station
t2 = 120.0  # Time to reach second station
t3 = 180.0  # Time to reach final station

# Desired stopping times at intermediate stations (passenger boarding/alighting)
stop_time_1 = 15.0  # seconds at station 1
stop_time_2 = 15.0  # seconds at station 2

# Time transformation function 
function phi(s, t_start, t_end)
    return t_start + (t_end - t_start) * s
end
nothing # hide

```
Model defintion 
```@example main-cit

# Tram dynamics function
function tram_dynamics(t, position, velocity, control)
    # Simple point mass dynamics: F = ma
    # position' = velocity
    # velocity' = acceleration (control input)
    
    acceleration = control
    
    # Return [position', velocity']
    return [velocity, acceleration]
end
nothing # hide

```
```@example main-cit

@def ocp_3arc begin
    
    s ∈ [0, 1], time
    y ∈ R^6,    state    # [pos1, vel1, pos2, vel2, pos3, vel3] for each arc
    w ∈ R^3,    control  # [acc1, acc2, acc3] for each arc
    
    # Arc 1: Station 0 to Station 1
    x1 = y[1:2]
    pos1 = y[1]
    vel1 = y[2]
    acc1 = w[1]
    
    # Arc 2: Station 1 to Station 2 (with stop time)
    x2 = y[3:4]
    pos2 = y[3]
    vel2 = y[4]
    acc2 = w[2]
    
    # Arc 3: Station 2 to Station 3 (with stop time)
    x3 = y[5:6]
    pos3 = y[5]
    vel3 = y[6]
    acc3 = w[3]
    
    # State constraints for Arc 1
    0 ≤ pos1(s) ≤ 3000
    0 ≤ vel1(s) ≤ max_velocity
    max_deceleration ≤ acc1(s) ≤ max_acceleration
    
    # State constraints for Arc 2
    0 ≤ pos2(s) ≤ 3000
    0 ≤ vel2(s) ≤ max_velocity
    max_deceleration ≤ acc2(s) ≤ max_acceleration
    
    # State constraints for Arc 3
    0 ≤ pos3(s) ≤ 3000
    0 ≤ vel3(s) ≤ max_velocity
    max_deceleration ≤ acc3(s) ≤ max_acceleration
    
    # Initial conditions: Arc 1 (start from Station 0)
    pos1(0) == station_0
    vel1(0) == 0.0
    
    # Final conditions: Arc 1 (arrive at Station 1)
    pos1(1) == station_1
    vel1(1) == 0.0  # Must stop at station
    
    x2(0) - x1(1) == [0.0, 0.0]  # Position and velocity continuity
    
    # Final conditions: Arc 2 (arrive at Station 2)
    pos2(1) == station_2
    vel2(1) == 0.0  # Must stop at station
    
    x3(0) - x2(1) == [0.0, 0.0]  # Position and velocity continuity
    
    # Final conditions: Arc 3 (arrive at Station 3)
    pos3(1) == station_3
    vel3(1) == 0.0  # Stop at final destination
    
    # Dynamics for each arc
    ẏ(s) == [
        # Arc 1: t0 to t1
        (t1 - t0) * tram_dynamics(phi(s, t0, t1), pos1(s), vel1(s), acc1(s))...,
        # Arc 2: t1+stop_time_1 to t2  
        (t2 - t1 - stop_time_1) * tram_dynamics(phi(s, t1 + stop_time_1, t2), pos2(s), vel2(s), acc2(s))...,
        # Arc 3: t2+stop_time_2 to t3
        (t3 - t2 - stop_time_2) * tram_dynamics(phi(s, t2 + stop_time_2, t3), pos3(s), vel3(s), acc3(s))...
    ]
    
    # Objective: Minimize total energy consumption across all arcs
    ∫(acc1(s)^2 + acc2(s)^2 + acc3(s)^2) → min
    
end
nothing # hide

```
```@example main-cit

# Solve the three-arc problem
println("Solving three-arc tram problem...")
sol_3arc = solve(ocp_3arc; grid_size = 150, max_iter = 2000, tol = 1e-6, acceptable_tol = 1e-5)
println("Three-arc problem solved successfully!")
```
