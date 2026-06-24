# Flapping Wing 2D CFD Simulation

[![Julia](https://img.shields.io/badge/Julia-1.10-9558B2?logo=julia)](https://julialang.org)
[![WaterLily](https://img.shields.io/badge/WaterLily-1.7-007ACC)](https://github.com/weymouth/WaterLily.jl)

A 2D flapping wing simulation using the **Boundary Data Immersion Method (BDIM)** in [WaterLily.jl](https://github.com/weymouth/WaterLily.jl). The wing undergoes combined heave (vertical translation) and pitch (rotation) motion, representative of biological flight and bio-inspired propulsion.

<p align="center">
  <img src="docs/img/vorticity_placeholder.png" alt="Vorticity field visualization (add your ParaView screenshot here)" width="600"/>
  <br><em>Example vorticity field — replace with your own ParaView screenshot.</em>
</p>

---

## Physics

### Wing Kinematics

The foil center follows a combined heave–pitch trajectory with a $\pi/2$ phase offset:

$$
\begin{aligned}
h(t) &= h_{\text{amp}} \cdot c \cdot \sin(2\pi f t) \\
\theta(t) &= \theta_{\text{amp}} \cdot \cos(2\pi f t)
\end{aligned}
$$

where $c$ is the chord length, $f$ is the flapping frequency, $h_{\text{amp}}$ is the heave amplitude (in chords), and $\theta_{\text{amp}}$ is the pitch amplitude (in radians). Heave leads pitch by $90^\circ$, producing the characteristic flapping stroke.

### Dimensionless Numbers

**Reynolds number** — ratio of inertial to viscous forces:

$$
Re = \frac{U c}{\nu}
$$

**Strouhal number** — ratio of unsteady to inertial forces, based on total heave amplitude $A = 2 h_{\text{amp}} c$:

$$
St = \frac{f A}{U}
$$

The flapping frequency is derived from the Strouhal definition:

$$
f = \frac{St \cdot U}{2 h_{\text{amp}} c}
$$

### Force Coefficients

Drag and lift coefficients are computed at each timestep using WaterLily's built-in force integration (`WaterLily.total_force(sim)`). The coefficients are normalized as:

$$
C_d = \frac{F_x}{\frac{1}{2} \rho U^2 c}, \qquad
C_l = \frac{F_y}{\frac{1}{2} \rho U^2 c}
$$

where $\rho = 1$ in lattice units.

### Boundary Data Immersion Method

WaterLily uses the BDIM, which represents solid boundaries through a smooth kernel convolution rather than body-conformal meshes. The body is defined by:
- A **signed distance function (SDF)** $\phi(\mathbf{x}, t)$ — negative inside the body, positive outside
- A **map function** $\mathbf{\chi}(\mathbf{x}, t)$ — transforms world coordinates to body coordinates

This simulation uses a rounded-rectangle foil (NACA-like profile) with thickness $0.12c$.

### Smooth Startup

An impulsive start can generate transient vortices that contaminate early-time statistics. A smooth ramp-up envelope is applied over the first flapping period $T$:

$$
\text{ramp}(t) = \begin{cases}
\frac{1}{2}\left(1 - \cos\frac{\pi t}{T}\right), & t < T \\
1, & t \geq T
\end{cases}
$$

The heave and pitch amplitudes are multiplied by this envelope:
$h(t) = \text{ramp}(t) \cdot h_0(t)$, $\theta(t) = \text{ramp}(t) \cdot \theta_0(t)$.

> **Note on domain size:** The domain extends 9 chord lengths downstream of the foil, which is adequate for $Re = 250$ at moderate Strouhal numbers. For higher Reynolds numbers or Strouhal numbers, a larger wake region is recommended to prevent vortex reflections from the outlet boundary.

---

## Getting Started

### Prerequisites

- [Julia](https://julialang.org/downloads/) ≥ 1.10
- [ParaView](https://www.paraview.org/download/) (for visualization)

### Run the Simulation

```bash
git clone https://github.com/aryandeshmukh1510/FlappingWingSimulation
cd FlappingWingTest1

# Instantiate project dependencies
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Run with default parameters
julia --project=. src/FlappingWing.jl
```

### Custom Parameters

Default parameters can be overridden by either:

1. **Modifying the struct defaults** in `src/FlappingWing.jl`, or
2. **Creating a custom script** (outside the repo):

```julia
include("src/FlappingWing.jl")
using .FlappingWing

params = FlappingWingParams(
    chord = 64,       # higher resolution
    Re = 1000,        # higher Reynolds number
    St = 0.6,         # higher Strouhal number
    t_end = 20.0,     # longer simulation
)

run!(params; output_dir="output")
```

---

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `chord`   | 32      | Chord length in lattice units |
| `thk_ratio` | 0.12 | Foil thickness-to-chord ratio |
| `Re`      | 250     | Reynolds number |
| `U`       | 1.0     | Free-stream velocity (lu) |
| `St`      | 0.3     | Strouhal number |
| `h_amp`   | 1.0     | Heave amplitude (fraction of chord) |
| `θ_amp`   | π/4     | Pitch amplitude (radians) |
| `nx`      | 12      | Domain size in chords (x-direction) |
| `ny`      | 8       | Domain size in chords (y-direction) |
| `center_frac_x` | 0.25 | Foil pivot x-position (fraction of nx) |
| `center_frac_y` | 0.50 | Foil pivot y-position (fraction of ny) |
| `t_end`   | 10.0    | Simulation end time (convective units) |
| `write_interval` | 0.1 | Time between VTK exports |
| `Δt_max`  | 0.25    | Max time step (CFL constraint) |

---

## Visualization in ParaView

1. Open ParaView
2. **File → Open** → select `output/flapping_sim.pvd`
3. Click **Apply** in the Properties panel
4. In the dropdown next to the toolbar, select a field to visualize:
   - `u` — velocity magnitude (good overview)
   - `p` — pressure field (shows high/low regions)
   - `body_sdf` — signed distance function (shows the body)
5. To view **vorticity**:
   - Select the pipeline, go to **Filters → Alphabetical → Gradient Of Unstructured Dataset**
   - Set **Scalar Array** to `u` and **Result Arrays** to `Vorticity`
   - Click **Apply**
6. Use **Play** (spacebar) to animate through timesteps
7. For publication-quality output:
   - Use the **Save Animation** or **Export Scene** tools
   - Recommended: vorticity colormap with a diverging palette ("Cool to Warm")

---

## Project Structure

```
FlappingWingTest1/
├── src/
│   └── FlappingWing.jl      # Simulation module (struct, SDF, motion, force output)
├── output/                   # VTK snapshots (gitignored)
│   ├── flapping_sim.pvd      # Time series file for ParaView
│   └── flapping_sim_*.vti    # Individual timestep files
├── docs/
│   └── img/                  # Screenshots and figures
├── .gitignore
├── LICENSE                   # MIT License
├── Project.toml              # Julia project dependencies
└── README.md
```

---

## Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| [WaterLily.jl](https://github.com/weymouth/WaterLily.jl) | 1.7 | Immersed boundary CFD solver |
| [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) | 1.9 | Stack-allocated vector math |
| [WriteVTK.jl](https://github.com/jipolanco/WriteVTK.jl) | 1.21 | VTK file export for ParaView |
| [Printf.jl](https://docs.julialang.org/en/v1/stdlib/Printf/) | (stdlib) | Formatted console output |

---

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.

Copyright (c) 2026 Aryan Deshmukh
