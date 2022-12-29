export AbstractXie, Xie, XieNoisy, xie_affect!, xie_turnrate

abstract type AbstractXie{D} <: AbstractMicrobe{D} end

Base.@kwdef mutable struct Xie{D} <: AbstractXie{D}
    id::Int
    pos::NTuple{D,Float64} = ntuple(zero, D)
    motility = RunReverseFlick(speed_forward = Degenerate(46.5))
    vel::NTuple{D,Float64} = rand_vel(D, motility)
    turn_rate_forward::Float64 = 2.3 # 1/s
    turn_rate_backward::Float64 = 1.9 # 1/s
    state_m::Float64 = 0.0 # s
    state_z::Float64 = 0.0 # s
    state::Float64 = 0.0 # s
    adaptation_time_m::Float64 = 1.29 # s
    adaptation_time_z::Float64 = 0.28 # s 
    gain_forward::Float64 = 2.7 # 1/s
    gain_backward::Float64 = 1.6 # 1/s
    binding_affinity::Float64 = 0.39 # μM
    rotational_diffusivity::Float64 = 0.26 # rad²/s
    radius::Float64 = 0.5 # μm
end

Base.@kwdef mutable struct XieNoisy{D} <: AbstractXie{D}
    id::Int
    pos::NTuple{D,Float64} = ntuple(zero, D)
    motility = RunReverseFlick(speed_forward = Degenerate(46.5))
    vel::NTuple{D,Float64} = rand_vel(D, motility)
    turn_rate_forward::Float64 = 2.3 # 1/s
    turn_rate_backward::Float64 = 1.9 # 1/s
    state_m::Float64 = 0.0 # s
    state_z::Float64 = 0.0 # s
    state::Float64 = 0.0 # s
    adaptation_time_m::Float64 = 1.29 # s
    adaptation_time_z::Float64 = 0.28 # s 
    gain_forward::Float64 = 2.7 # 1/s
    gain_backward::Float64 = 1.6 # 1/s
    binding_affinity::Float64 = 0.39 # μM
    chemotactic_precision::Float64 = 6.0 # 1
    rotational_diffusivity::Float64 = 0.26 # rad²/s
    radius::Float64 = 0.5 # μm
end

function xie_affect!(microbe::Xie, model)
    Δt = model.timestep
    c = model.concentration_field(microbe.pos, model)
    K = microbe.binding_affinity
    ϕ = log(1.0 + c/K)
    τ_m = microbe.adaptation_time_m
    τ_z = microbe.adaptation_time_z
    a₀ = (τ_m*τ_z)/(τ_m - τ_z)
    m = microbe.state_m
    z = microbe.state_z
    m += (ϕ - m/τ_m)*Δt
    z += (ϕ - z/τ_z)*Δt
    microbe.state_m = m
    microbe.state_z = z
    microbe.state = a₀ * (m/τ_m - z/τ_z)
    return nothing
end

function xie_affect!(microbe::XieNoisy, model; ε=1e-16)
    Δt = model.timestep
    Dc = model.compound_diffusivity
    c = model.concentration_field(microbe.pos, model)
    K = microbe.binding_affinity
    a = microbe.radius
    Π = microbe.chemotactic_precision
    σ = Π * 0.04075 * sqrt(3*c / (5*π*Dc*a*Δt))
    M = rand(Normal(c,σ))
    ϕ = log(1.0 + max(M/K, -1+ε))
    τ_m = microbe.adaptation_time_m
    τ_z = microbe.adaptation_time_z
    a₀ = (τ_m*τ_z)/(τ_m - τ_z)
    m = microbe.state_m
    z = microbe.state_z
    m += (ϕ - m/τ_m)*Δt
    z += (ϕ - z/τ_z)*Δt
    microbe.state_m = m
    microbe.state_z = z
    microbe.state = a₀ * (m/τ_m - z/τ_z)
    return nothing
end

function xie_turnrate(microbe, model)
    if microbe.motility isa AbstractMotilityTwoStep
        return xie_turnrate_twostep(microbe, model)
    else
        return xie_turnrate_onestep(microbe, model)
    end
end
function xie_turnrate_onestep(microbe, model)
    S = microbe.state
    ν₀ = microbe.turn_rate_forward
    β = microbe.gain_forward
    return ν₀*(1 + β*S)
end
function xie_turnrate_twostep(microbe, model)
    S = microbe.state
    if motility.state == Forward
        ν₀ = microbe.turn_rate_forward
        β = microbe.gain_forward
    elseif motility.state == Backward
        ν₀ = microbe.turn_rate_backward
        β = microbe.gain_backward
    end
    return ν₀*(1 + β*S)
end

function microbe_step!(microbe::AbstractXie, model::ABM)
    microbe_step!(
        microbe, model;
        affect! = xie_affect!,
        turnrate = xie_turnrate
    )
end