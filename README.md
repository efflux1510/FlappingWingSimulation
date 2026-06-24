# Flapping Wing 2D CFD Simulation

[![Julia](https://img.shields.io/badge/Julia-1.10-9558B2?logo=julia)](https://julialang.org)
[![WaterLily](https://img.shields.io/badge/WaterLily-1.7-007ACC)](https://github.com/weymouth/WaterLily.jl)

I built this to simulate a 2D flapping wing — like a bird or fish — using [WaterLily.jl](https://github.com/weymouth/WaterLily.jl), a CFD solver that handles complex shapes without needing a body-fitted mesh. The wing does a combined up-down (heave) and tilt (pitch) motion. It prints force coefficients (Cd, Cl) as it runs, and exports everything to ParaView so you can watch the vortices form and shed.

---

## Physics

### Wing motion

The foil moves like this:

$$
\begin{aligned}
h(t) &= h_{\text{amp}} \cdot c \cdot \sin(2\pi f t) \\
\theta(t) &= \theta_{\text{amp}} \cdot \cos(2\pi f t)
\end{aligned}
$$

Heave ($h$) and pitch ($\theta$) are $90^\circ$ out of phase — the wing plunges up and down while tilting, just like a real flapping wing. $c$ is the chord length, $f$ is the flapping frequency.

### Reynolds and Strouhal numbers

Two numbers control the flow physics:

$$
Re = \frac{U c}{\nu}, \qquad
St = \frac{f A}{U}
$$

$Re$ compares inertia to viscosity (laminar vs turbulent), $St$ compares flapping speed to flow speed. $A = 2 h_{\text{amp}} c$ is the total heave amplitude. The frequency comes from the Strouhal definition:

$$
f = \frac{St \cdot U}{2 h_{\text{amp}} c}
$$

### Force coefficients

At every timestep, WaterLily integrates the pressure and viscous forces on the body and I normalize them:

$$
C_d = \frac{F_x}{\frac{1}{2} \rho U^2 c}, \qquad
C_l = \frac{F_y}{\frac{1}{2} \rho U^2 c}
$$

($\rho = 1$ in lattice units.) You'll see these printed in the console as the simulation runs.

### How the body works

WaterLily uses the **Boundary Data Immersion Method (BDIM)** — no mesh around the body, just a signed distance function (SDF) that says "inside" vs "outside", and a map function that moves coordinates from the world frame to the body frame. The foil here is a rounded rectangle, basically a flat plate with rounded edges, 12% thick relative to chord.

### Smooth startup

If you just start flapping from rest, you get a big lurch in the flow that takes a while to wash out. I added a smooth ramp over the first period:

$$
\text{ramp}(t) = \begin{cases}
\frac{1}{2}\left(1 - \cos\frac{\pi t}{T}\right), & t < T \\
1, & t \geq T
\end{cases}
$$

The heave and pitch get multiplied by this, so they fade in gently.

> **One caveat:** The domain is 9 chord lengths behind the foil. That's fine for $Re = 250$ at moderate Strouhal numbers, but if you crank up $Re$ or $St$, vortices might bounce off the outlet boundary. Just make the domain bigger.

---

## Getting Started

### You'll need

- [Julia](https://julialang.org/downloads/) ≥ 1.10
- [ParaView](https://www.paraview.org/download/) to look at the results

### Run it

```bash
git clone https://github.com/aryandeshmukh1510/FlappingWingSimulation
cd FlappingWingTest1

# Download dependencies
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Run with defaults
julia --project=. src/FlappingWing.jl
```

### Change parameters

Tweak anything by passing keyword arguments:

```julia
include("src/FlappingWing.jl")
using .FlappingWing

params = FlappingWingParams(
    chord = 64,       # higher resolution
    Re = 1000,        # more interesting flow
    St = 0.6,         # faster flapping
    t_end = 20.0,     # longer sim
)

run!(params; output_dir="output")
```

Or just edit the defaults in `src/FlappingWing.jl` directly.

---

## Parameters

All the knobs and what they do:

| Parameter | Default | What it is |
|-----------|---------|------------|
| `chord`   | 32      | Chord length in grid cells |
| `thk_ratio` | 0.12 | Foil thickness ÷ chord |
| `Re`      | 250     | Reynolds number |
| `U`       | 1.0     | Flow speed |
| `St`      | 0.3     | Strouhal number |
| `h_amp`   | 1.0     | Heave amplitude (× chord) |
| `θ_amp`   | π/4     | Pitch amplitude (radians) |
| `nx`      | 12      | Domain width in chords |
| `ny`      | 8       | Domain height in chords |
| `center_frac_x` | 0.25 | Foil x-position (fraction of width) |
| `center_frac_y` | 0.50 | Foil y-position (fraction of height) |
| `t_end`   | 10.0    | How long to simulate |
| `write_interval` | 0.1 | Time between saved frames |
| `Δt_max`  | 0.25    | Max timestep size |

---

## Viewing in ParaView

1. Open ParaView
2. **File → Open** → pick `output/flapping_sim.pvd`
3. Hit **Apply**
4. Pick a field to color by: `u` (velocity), `p` (pressure), or `body_sdf`
5. Want **vorticity**? **Filters → Alphabetical → Gradient Of Unstructured Dataset**, set **Scalar Array** to `u`, hit **Apply**
6. Hit **Play** (spacebar) to watch it run through time

---

## Project layout

```
FlappingWingTest1/
├── src/
│   └── FlappingWing.jl      # The main code
├── output/                   # VTK files (gitignored)
│   ├── flapping_sim.pvd
│   └── flapping_sim_*.vti
├── docs/
│   └── img/                  # For your ParaView screenshots
├── .gitignore
├── LICENSE
├── Project.toml
└── README.md
```

---

## Dependencies

| Package | What it does |
|---------|--------------|
| [WaterLily.jl](https://github.com/weymouth/WaterLily.jl) | The CFD engine |
| [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) | Fast fixed-size arrays |
| [WriteVTK.jl](https://github.com/jipolanco/WriteVTK.jl) | ParaView file export |
| Printf | Console formatting (built into Julia) |

---

## License

MIT — see [LICENSE](LICENSE).

Copyright (c) 2026 Aryan Deshmukh
