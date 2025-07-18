# [constraints-at-intermediate-times](@id tutorial-for-constraints-at-intermediate-times)
In this tutorial, we will explore an example for constraints at intermediate times
We want to showcase use cases where we have some intermediate times with constraints associated to each one:


## first simple example to modelize the orbit of the moon around the earth
Necessary imports  : 

```@example main-cit
#import Pkg; Pkg.add("SPICE")
#import Pkg; Pkg.add("DataInterpolations")
using SPICE
using OptimalControl
using NLPModelsIpopt
using DataInterpolations
using Downloads: download
using Plots
using Plots.PlotMeasures # for leftmargin, bottommargin
#using CSV
#using JLD2

# SPICE files that you have to download but only once
const LSK = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls"
const SPK = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp"
const PCK = "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc"

# Download kernels
download(LSK, "naif0012.tls")
download(SPK, "de440.bsp")
download(PCK, "pck00010.tpc")

furnsh("naif0012.tls") # Leap seconds kernel
furnsh("de440.bsp") # Ephemeris kernel
furnsh("pck00010.tpc")  # PCK kernel for planetary constants
```
```@example main-cit

# Astrodynamical parameters
G = 6.6743e-11 / 1000^3 # km^3 kg^-1 s^-2
M_moon = 7.346e22 # kg
M_earth = 5.972168e24 # kg
M_sun = 1.989e30 # kg
AU = 1.495978707e11 / 1e3 # km
LU = 384399 # km
TU = sqrt(LU^3 / G / (M_moon + M_earth))
MUnit =  M_moon + M_earth
mu_Earth = 3.9860044188e14 / 1000^3 / LU^3 * TU^2
mu_Moon = 4.90486959e12 / 1000^3 / LU^3 * TU^2
mu_Sun = 1.327124400189e20  / 1000^3 / LU^3 * TU^2
mass_moon = M_moon / MUnit
mass_earth = M_earth / MUnit
mass_sun = M_sun / MUnit

# Spacecraft parameters
g0 = 9.80665 / 1000 / LU * TU^2
Isp = 1660 / TU 
wet_mass = 344 / MUnit 
fuel_mass = 90 / MUnit 
dry_mass = (wet_mass - fuel_mass) 
wet_by_dry = (dry_mass + fuel_mass) / dry_mass
Tmax =  0.090 / 1000 / MUnit / LU * TU^2 #thrust force in N
Q = Tmax / Isp / g0

mass_scaling_factor = wet_mass

# Mission data that I just copy paste here for simplicity
t0_J200 = 1.107747244e9
tf_J200 = 1.112803048e9
N = 15000 
t_values = LinRange(t0_J200, tf_J200, N-3)
t_values = vcat(-5e9, 0.0, collect(t_values), tf_J200 * 5)
states_moon = zeros(6, N)
states_earth = zeros(6, N)
states_sun = zeros(6, N)
t0 = t0_J200 / TU
tf = tf_J200 / TU

# Interpolations of the bodies positions/velocities
for i = 1:N
    states_moon[:,i], _ = spkezr("MOON", t_values[i], "ECLIPJ2000", "NONE", "EARTH MOON BARYCENTER")
    states_earth[:,i], _ = spkezr("EARTH", t_values[i], "ECLIPJ2000", "NONE", "EARTH MOON BARYCENTER")
    states_sun[:,i], _ = spkezr("SUN", t_values[i], "ECLIPJ2000", "NONE", "EARTH MOON BARYCENTER")
end

# Interpolation of the position and velocity of the celestial bodies (Sun, Moon, Earth)
x_moon = DataInterpolations.LinearInterpolation(states_moon[1,:] / LU, t_values / TU) 
y_moon = DataInterpolations.LinearInterpolation(states_moon[2,:] / LU, t_values / TU) 
z_moon = DataInterpolations.LinearInterpolation(states_moon[3,:] / LU, t_values / TU) 

x_earth = DataInterpolations.LinearInterpolation(states_earth[1,:] / LU, t_values / TU) 
y_earth = DataInterpolations.LinearInterpolation(states_earth[2,:] / LU, t_values / TU) 
z_earth = DataInterpolations.LinearInterpolation(states_earth[3,:] / LU, t_values / TU) 

x_sun = DataInterpolations.LinearInterpolation(states_sun[1,:] / LU, t_values / TU) 
y_sun = DataInterpolations.LinearInterpolation(states_sun[2,:] / LU, t_values / TU) 
z_sun = DataInterpolations.LinearInterpolation(states_sun[3,:] / LU, t_values / TU) 

# Times
t0 = 2952.5065962806448
t1 = 2953.955962786101
t2 = 2954.4423792964794
t3 = 2955.398985817685
t4 = 2957.4519685112937
t5 = 2959.0986010622196
t6 = 2960.6895789585565

# States at the previous times
ap0      =  [-0.992463207586249; 0.7937735622182638; 0.07902677416947182; 0.4497665242298916;  0.5768992856452368; -0.08990077540870026]
occult1  =  [0.024011456548615025; 1.0839463211562548; -0.06721110889958525; 0.8370532723530986; -0.2849379073937811; -0.09055054145691602]
lunar2   =  [0.42191676845234877; 0.8365166413216278; -0.10058044976319964; 1.0641134073373073;  -0.7420861511380829; 0.2221102307286784]
ap3      =  [0.8866648294530993; 0.07703956973358668; 0.28758014926093095; 0.010909379448810359; -1.0104415437867642; 0.3065975859809964]
lunar4   =  [-0.5240519542473246; -0.8727170323241018; 0.06983526263136841; -1.0656424593386877; 0.4013541226489389; -0.24177176981236737]
ap5      =  [-1.2845625006191574; 0.07452729181167739; -0.02061966816776342; 0.03379116945381267; 0.7236740339189934; -0.05964259399261467]
occult6  =  [-0.5063291057829847; 0.9375405520827109; -0.0856392904912281; 0.906825723975573; 0.17758145410610057; -0.009952663265456345]

# Changement de variable pour le temps
function phi(s, t0, t1)
    return t0 + (t1 - t0) * s
end

# Dynamique
function Fu(t, r, v, m, u)
    
    # Earth & Moon & Sun gravity
    m = m * mass_scaling_factor

    x_e = [x_earth(t), y_earth(t), z_earth(t)]
    x_m = [x_moon(t), y_moon(t), z_moon(t)]
    x_s = [x_sun(t), y_sun(t), z_sun(t)]

    normrE = sqrt(sum((x_e - r).^2))
    normrM = sqrt(sum((x_m - r).^2))
    normrS = sqrt(sum((x_s - r).^2))

    normrSE = sqrt(sum((x_e - x_s).^2))
    normrMS = sqrt(sum((x_m - x_s).^2))

    dv = mu_Earth / normrE^3 * (x_e - r) + 
         mu_Moon  / normrM^3 * (x_m - r) + 
         mu_Sun   / normrS^3 * (x_s - r) - 
         (1 / (mass_earth + mass_moon)) * ( 
            x_s * (mu_Sun  * mass_earth / normrSE^3 + mu_Moon * mass_sun / normrMS^3) - 
            x_m * (mu_Moon * mass_sun   / normrMS^3) - 
            x_e * (mu_Sun * mass_earth / normrSE^3)
            )
    dvu = Tmax * u / m
    dv = dv + dvu

    dm = - Q * sqrt(sum(u.^2)) / mass_scaling_factor

    dx = [v..., dv..., dm]
    return dx

end
```


```@example main-cit
# Problem with 1 arc

@def ocp1arc begin

    s ∈ [0, 1], time
    y ∈ R^7,    state
    w ∈ R^3,    control

    r1 = y[1:3]
    v1 = y[4:6]
    m1 = y[7]
    u1 = w[1:3]

    [-5, -5, -5] ≤ r1(s) ≤ [5, 5, 5]
    [-5, -5, -5] ≤ v1(s) ≤ [5, 5, 5]
    0.5 ≤ m1(s) ≤ 1 # careful here ! dry mass > 0

    0 ≤ sum(u1(s).^2) ≤ 1, eq_u1

    r1(0) == occult1[1:3]
    v1(0) == occult1[4:6]
    m1(0) == 1
    r1(1) == lunar2[1:3]
    v1(1) == lunar2[4:6]

    ẏ(s) == (t2 - t1) * Fu(phi(s, t1, t2), r1(s), v1(s), m1(s), u1(s))

    m1(1) → max

end

sol = solve(ocp1arc; grid_size = 100, max_iter = 3000, tol = 1e-6, acceptable_tol = 1e-5)
sol = solve(ocp1arc; grid_size = 500, max_iter = 3000, tol = 1e-6, acceptable_tol = 1e-5, init = sol)
sol = solve(ocp1arc; grid_size = 900, max_iter = 3000, tol = 1e-6, acceptable_tol = 1e-5, init = sol)
plot(sol; state_bounds_style=:none, leftmargin=10mm, bottommargin=5mm)
```

## Example to modelize a trein passing through two stations (3 arcs)


