export Celani, CelaniNoisy, celani_affect!, celani_turnrate, microbe_step!

abstract type AbstractCelani{D} <: AbstractMicrobe{D} end

"""
    Celani{D} <: <: AbstractCelani{D} <: AbstractMicrobe{D}
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
Base.@kwdef mutable struct Celani{D} <: AbstractCelani{D}
    id::Int
    pos::NTuple{D,Float64} = ntuple(zero, D)
    motility = RunTumble(speed = Degenerate(30.0))
    vel::NTuple{D,Float64} = rand_vel(D) .* rand(motility.speed) # μm/s
    turn_rate::Float64 = 1/0.67 # 1/s
    state::Vector{Float64} = [0.,0.,0.,1.] # 1
    rotational_diffusivity = 0.26 # rad²/s
    gain::Float64 = 50.0 # 1
    memory::Float64 = 1.0 # s
    radius::Float64 = 0.5 # μm
end # struct

"""
    CelaniNoisy{D} <: AbstractCelani{D} <: AbstractMicrobe{D}
Model of chemotactic bacterium using the response kernel from 'Celani and Vergassola (2010) PNAS'
with a custom implementation of noisy sensing based on 'Berg and Purcell (1977) Biophys J'.
The intensity of noise implicit in bacterial sensing is represented by a dimensionless
factor `chemotactic_precision` (same approach used by 'Brumley et al. (2019) PNAS').
`chemotactic_precision = 1` represents the minimum theoretical noise.

Default parameters:
- motility = RunTumble(speed = Degenerate(30.0))
- turn_rate = 1.49 Hz
- state = zeros(4)
- rotational_diffusivity = 0.26 rad²/s
- gain = 50.0
- memory = 1 s
- chemotactic_precision = 1.0
- radius = 0 μm
"""
Base.@kwdef mutable struct CelaniNoisy{D} <: AbstractCelani{D}
    id::Int
    pos::NTuple{D,Float64} = ntuple(zero, D)
    motility = RunTumble(speed = Degenerate(30.0))
    vel::NTuple{D,Float64} = rand_vel(D) .* rand(motility.speed) # μm/s
    turn_rate::Float64 = 1/0.67 # 1/s
    state::Vector{Float64} = [0.,0.,0.,1.] # 1
    rotational_diffusivity = 0.26 # rad²/s
    gain::Float64 = 50.0 # 1
    memory::Float64 = 1.0 # s
    chemotactic_precision::Float64 = 1.0 # 1
    radius::Float64 = 0.5 # μm
end # struct

function celani_affect!(microbe::Celani, model)
    Δt = model.timestep
    u = model.concentration_field(microbe.pos, model)
    γ = microbe.memory
    λ = 1/γ
    β = microbe.gain
    S = microbe.state
    S[1] = S[1] + (-λ*S[1] + u)*Δt
    S[2] = S[2] + (-λ*S[2] + S[1])*Δt
    S[3] = S[3] + (-λ*S[3] + 2*S[2])*Δt
    S[4] = 1 - β*(λ^2*S[2] - λ^3/2*S[3])
    return nothing
end # function

function celani_affect!(microbe::CelaniNoisy, model)
    Δt = model.timestep
    Dc = model.compound_diffusivity
    u = model.concentration_field(microbe.pos, model)
    a = microbe.radius
    Π = microbe.chemotactic_precision
    σ = Π * 0.04075 * sqrt(3*u / (5*π*Dc*a*Δt)) # noise (Berg-Purcell)
    M = rand(Normal(u,σ)) # measurement
    γ = microbe.memory
    λ = 1/γ
    β = microbe.gain
    S = microbe.state
    S[1] = S[1] + (-λ*S[1] + M)*Δt
    S[2] = S[2] + (-λ*S[2] + S[1])*Δt
    S[3] = S[3] + (-λ*S[3] + 2*S[2])*Δt
    S[4] = 1 - β*(λ^2*S[2] - λ^3/2*S[3])
    return nothing
end # function

function celani_turnrate(microbe, model)
    ν₀ = microbe.turn_rate # unbiased
    S = microbe.state[4]
    return ν₀*S # modulated turn rate
end # function

function microbe_step!(microbe::AbstractCelani, model)
    microbe_step!(
        microbe, model;
        affect! = celani_affect!,
        turnrate = celani_turnrate
    )
end # function