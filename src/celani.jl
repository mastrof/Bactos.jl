export MicrobeCelani, celani_affect!, celani_turnrate, microbe_step!

"""
    MicrobeCelani{D} <: AbstractMicrobe{D}
Model of chemotactic bacterium using the response kernel from 'Celani and Vergassola (2010) PNAS'.
Default parameters:
- motility = RunTumble(speed = Degenerate(30.0))
- turn_rate = 1.49 Hz
- state = zeros(4)
- rotational_diffusivity = 0.26 rad²/s
- gain = 50.0
- memory = 1 s
- radius = 0 μm
"""
Base.@kwdef mutable struct MicrobeCelani{D} <: AbstractMicrobe{D}
    id::Int
    pos::NTuple{D,Float64} = ntuple(zero, D)
    motility = RunTumble(speed = Degenerate(30.0))
    vel::NTuple{D,Float64} = rand_vel(D) .* rand(motility.speed) # μm/s
    turn_rate::Float64 = 1/0.67 # 1/s
    state::Vector{Float64} = zeros(4) # 1
    rotational_diffusivity = 0.26 # rad²/s
    gain::Float64 = 50.0 # 1
    memory::Float64 = 1.0 # s
    radius::Float64 = 0.0 # μm
end # struct

function celani_affect!(microbe, model)
    Δt = model.timestep
    u = model.concentration_field(microbe.pos, model)
    ∇u = model.concentration_gradient(microbe.pos, model)
    ∂ₜu = model.concentration_time_derivative(microbe.pos, model)
    du_dt = dot(microbe.vel, ∇u) + ∂ₜu
    γ = microbe.memory
    _γ = 1/γ
    β = microbe.gain
    S = microbe.state
    S[1] = Δt * (-S[1]*_γ + du_dt)
    S[2] = Δt * (S[1]*_γ - S[2]*_γ)
    S[3] = Δt * (S[2]*2*_γ - S[3]*_γ)
    S[4] = (1 - β*S[3])
    return nothing
end # function

function celani_turnrate(microbe, model)
    ν₀ = microbe.turn_rate # unbiased
    S = microbe.state[4]
    return ν₀*S # modulated turn rate
end # function

function microbe_step!(microbe::MicrobeCelani, model)
    microbe_step!(
        microbe, model;
        affect! = celani_affect!,
        turnrate = celani_turnrate
    )
end # function