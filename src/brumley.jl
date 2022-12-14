export MicrobeBrumley, brumley_affect!, brumley_turnrate


"""
    MicrobeBrumley{D} <: AbstractMicrobe{D}
Model of chemotactic bacterium from 'Brumley et al. (2019) PNAS'.
The model is optimized for simulation of marine bacteria and accounts
for the presence of (gaussian) sensing noise in the chemotactic pathway.
Default parameters:
- motility = RunReverseFlick(speed_forward = Degenerate(46.5))
- turn_rate = 2.22 Hz → '1/τ₀'
- state = 0.0 → 'S'
- rotational_diffusivity = 0.035 rad²/S
- adaptation_time = 1.3 s → 'τₘ'
- receptor_gain = 50.0 μM⁻¹ → 'κ'
- motor_gain = 50.0 → 'Γ'
- chemotactic_precision = 6.0 → 'Π'
- radius = 0.5 μm → 'a'
"""
Base.@kwdef mutable struct MicrobeBrumley{D} <: AbstractMicrobe{D}
    id::Int
    pos::NTuple{D,Float64} = ntuple(zero, D)
    motility = RunReverseFlick(speed_forward = Degenerate(46.5))
    vel::NTuple{D,Float64} = rand_vel(D, motility) # μm/s
    turn_rate::Float64 = 1/0.45 # 1/s
    state::Float64 = 0.0
    rotational_diffusivity::Float64 = 0.035 # rad²/s
    adaptation_time::Float64 = 1.3 # s
    receptor_gain::Float64 = 50.0 # 1/μM
    motor_gain::Float64 = 50.0 # 1
    chemotactic_precision::Float64 = 6.0 # 1
    radius::Float64 = 0.5 # μm
end # struct

function brumley_affect!(microbe, model)
    Δt = model.timestep
    Dc = model.compound_diffusivity
    τₘ = microbe.adaptation_time
    α = exp(-Δt/τₘ) # memory persistence factor
    a = microbe.radius
    Π = microbe.chemotactic_precision
    κ = microbe.receptor_gain
    u = model.concentration_field(microbe.pos, model)
    ∇u = model.concentration_gradient(microbe.pos, model)
    ∂ₜu = model.concentration_time_derivative(microbe.pos, model)
    # gradient measurement
    μ = dot(microbe.vel, ∇u) + ∂ₜu # mean
    σ = Π * 0.04075 * sqrt(3*u / (π*a*Dc*Δt^3)) # noise
    M = rand(Normal(μ,σ)) # measurement
    # update internal state
    S = microbe.state
    microbe.state = α*S + (1-α)*κ*τₘ*M
    return nothing
end # function

function brumley_turnrate(microbe, model)
    ν₀ = microbe.turn_rate # unbiased
    Γ = microbe.motor_gain
    S = microbe.state
    return (1 + exp(-Γ*S)) * ν₀/2 # modulated turn rate
end # function

function microbe_step!(microbe::MicrobeBrumley, model)
    microbe_step!(
        microbe, model;
        affect! = brumley_affect!,
        turnrate = brumley_turnrate
    )
end # function