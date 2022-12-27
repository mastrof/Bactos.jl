export MicrobeBrownBerg, brownberg_affect!, brownberg_turnrate, microbe_step!


"""
    MicrobeBrownBerg{D} <: AbstractMicrobe{D}
Model of chemotactic E.coli from 'Brown and Berg (1974) PNAS'.
Default parameters:
- motility = RunTumble(speed = Degenerate(30.0))
- turn_rate = 1.49 Hz
- state = 0.0 → corresponds to "weighted dPb/dt" in the paper
- rotational_diffusivity = 0.035 rad²/s
- motor_gain = 660 s
- receptor_binding_constant = 100 μM
- adaptation_time = 1 s
- radius = 0 μm
"""
Base.@kwdef mutable struct MicrobeBrownBerg{D} <: AbstractMicrobe{D}
    id::Int
    pos::NTuple{D,Float64} = ntuple(zero, D)
    motility = RunTumble(speed = Degenerate(30.0))
    vel::NTuple{D,Float64} = rand_vel(D, motility) # μm/s
    turn_rate::Float64 = 1/0.67 # 1/s
    state::Float64 = 0.0 # 1
    rotational_diffusivity = 0.035 # rad²/s
    motor_gain::Float64 = 660.0 # s
    receptor_binding_constant::Float64 = 100.0 # μM
    adaptation_time::Float64 = 1.0 # s
    radius::Float64 = 0.0 # μm
end # struct

function brownberg_affect!(microbe, model)
    Δt = model.timestep
    τₘ = microbe.adaptation_time
    β = Δt / τₘ # memory loss factor
    KD = microbe.receptor_binding_constant
    S = microbe.state # weighted dPb/dt at previous step
    u = model.concentration_field(microbe.pos, model)
    ∇u = model.concentration_gradient(microbe.pos, model)
    ∂ₜu = model.concentration_time_derivative(microbe.pos, model)
    du_dt = dot(microbe.vel, ∇u) + ∂ₜu
    M = KD / (KD + u)^2 * du_dt # dPb/dt from new measurement
    microbe.state = β*M + S*exp(-β) # new weighted dPb/dt
    return nothing
end # function

function brownberg_turnrate(microbe, model)
    ν₀ = microbe.turn_rate # unbiased
    g = microbe.motor_gain
    S = microbe.state
    return ν₀*exp(-g*S) # modulated turn rate
end # function

function microbe_step!(microbe::MicrobeBrownBerg, model::ABM)
    microbe_step!(
        microbe, model;
        affect! = brownberg_affect!,
        turnrate = brownberg_turnrate
    )
end # function