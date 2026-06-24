module FlappingWing

using WaterLily
using StaticArrays
using WriteVTK
using Printf

export FlappingWingParams, run!

Base.@kwdef mutable struct FlappingWingParams
    chord::Int = 32
    thk_ratio::Float64 = 0.12
    Re::Float64 = 250.0
    U::Float64 = 1.0
    St::Float64 = 0.3
    h_amp::Float64 = 1.0
    θ_amp::Float64 = π / 4
    nx::Int = 12
    ny::Int = 8
    center_frac_x::Float64 = 0.25
    center_frac_y::Float64 = 0.5
    t_end::Float64 = 10.0
    write_interval::Float64 = 0.1
    Δt_max::Float64 = 0.25
end

_period(p) = 2 * p.h_amp * p.chord / (p.St * p.U)
_ang_frequency(p) = 2π / _period(p)
_center(p) = SA[p.center_frac_x * p.nx * p.chord, p.center_frac_y * p.ny * p.chord]
_viscosity(p) = p.U * p.chord / p.Re

function make_sdf(p::FlappingWingParams)
    thk = p.chord * p.thk_ratio
    return function sdf(x, t)
        y = x .- SA[clamp(x[1], -p.chord / 2, p.chord / 2), 0.0]
        return √sum(abs2, y) - thk / 2
    end
end

function make_motion(p::FlappingWingParams)
    ω = _ang_frequency(p)
    T = _period(p)
    center = _center(p)

    return function map(x, t)
        ramp = t < T ? 0.5 * (1 - cos(π * t / T)) : 1.0

        h = ramp * p.h_amp * p.chord * sin(ω * t)
        θ = ramp * p.θ_amp * cos(ω * t)

        R_inv = SA[cos(θ)  sin(θ); -sin(θ)  cos(θ)]

        x_shifted = x - center - SA[0.0, h]
        return R_inv * x_shifted
    end
end

function make_simulation(p::FlappingWingParams)
    body = AutoBody(make_sdf(p), make_motion(p))
    nx_lu, ny_lu = p.nx * p.chord, p.ny * p.chord
    return Simulation((nx_lu, ny_lu), (p.U, 0.0), p.chord; ν=_viscosity(p), body)
end

function run!(p::FlappingWingParams=FlappingWingParams(); output_dir="output")
    sim = make_simulation(p)
    ω = _ang_frequency(p)
    T = _period(p)
    center = _center(p)
    f = 1 / T
    n_periods = p.t_end / T

    mkpath(output_dir)

    writer = vtkWriter("flapping_sim"; dir=output_dir, attrib=Dict(
        "p" => sim -> sim.flow.p,
        "u" => sim -> sim.flow.u,
        "body_sdf" => sim -> sim.flow.σ,
    ))

    println("=" ^ 60)
    println("  Flapping Wing 2D Simulation")
    println("=" ^ 60)
    println("  Reynolds number:    Re = $(p.Re)")
    println("  Strouhal number:    St = $(p.St)")
    println("  Chord length:       c  = $(p.chord) lu")
    println("  Domain:             $(p.nx)×$(p.ny) chords ($(p.nx*p.chord)×$(p.ny*p.chord) lu)")
    println("  Foil center:        ($(round(center[1], digits=1)), $(round(center[2], digits=1))) lu")
    println("  Flapping frequency: f = $(round(f, digits=4))")
    println("  Flapping period:    T = $(round(T, digits=4))")
    println("  Sim duration:       t_end = $(p.t_end) ($(round(n_periods, digits=2)) periods)")
    println("  Viscosity:          ν = $(round(_viscosity(p), digits=6))")
    println("  Output directory:   $(abspath(output_dir))/")
    println("-" ^ 60)
    println("  Step   Time     Cd          Cl")
    println("-" ^ 60)

    forces = Float64[]
    step = 0

    for t in 0.0:p.write_interval:p.t_end
        sim_step!(sim, t; remeasure=true)
        save!(writer, sim)

        F = WaterLily.total_force(sim)
        Cd = F[1] / (0.5 * p.U^2 * p.chord)
        Cl = F[2] / (0.5 * p.U^2 * p.chord)

        push!(forces, Cd, Cl)
        step += 1

        println("  $(lpad(step, 3))   $(lpad(round(t, digits=2), 5))   $(@sprintf("%+9.5f", Cd))  $(@sprintf("%+9.5f", Cl))")
    end

    close(writer)

    n = length(forces) ÷ 2
    mean_Cd = sum(forces[1:2:end]) / n
    mean_Cl = sum(forces[2:2:end]) / n
    max_Cd = maximum(abs.(forces[1:2:end]))
    max_Cl = maximum(abs.(forces[2:2:end]))

    println("-" ^ 60)
    println("  Simulation complete.")
    println("  Steps written:      $(step)")
    println("  Periods simulated:  $(round(n_periods, digits=2))")
    println("  Mean Cd:            $(@sprintf("%+9.5f", mean_Cd))")
    println("  Mean Cl:            $(@sprintf("%+9.5f", mean_Cl))")
    println("  Max |Cd|:           $(round(max_Cd, digits=5))")
    println("  Max |Cl|:           $(round(max_Cl, digits=5))")
    println("  Output:             $(joinpath(output_dir, "flapping_sim.pvd"))")
    println("=" ^ 60)

    return (periods=n_periods, steps=step, mean_Cd=mean_Cd, mean_Cl=mean_Cl)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run!(FlappingWingParams())
end

end
