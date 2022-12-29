export
    AbstractMotility, AbstractMotilityOneStep, AbstractMotilityTwoStep,
    MotileState, TwoState, Forward, Backward, switch!,
    RunTumble, RunReverse, RunReverseFlick

"""
    AbstractMotility
General abstract interface for motility patterns.
"""
abstract type AbstractMotility end

"""
    AbstractMotilityOneStep
One-step motility patterns (`RunTumble`).
Subtypes must have at least the following fields:
- `speed`: distribution of microbe speed, new values extracted after each turn
- `polar`: distribution of polar angles
- `azimuthal`: distribution azimuthal angles
For 2-dimensional microbe types, only `polar` defines reorientations and `azimuthal` is ignored.
"""
abstract type AbstractMotilityOneStep <: AbstractMotility end

"""
    AbstractMotilityTwoStep
Two-step motility patterns (`RunReverse` and `RunReverseFlick`), with different
properties between forward and backward state of motion.
Subtypes must have at least the following fields:
- `speed_forward`: distribution of microbe speed, new values extracted after each turn
- `polar_forward`: distribution of in-plane reorientations for motile state 0
- `azimuthal_forward`: distribution of out-of-plane reorientations for motile state 0
- `speed_backward`: distribution of microbe speed, new values extracted after each turn
- `polar_backward`: distribution of in-plane reorientations for motile state 1
- `azimuthal_backward`: distribution of out-of-plane reorientations for motile state 1
- `motile_state`: defines current motile state (e.g. `Forward` or `Backward` for a `TwoState`)
For 2-dimensional microbe types, only `polar_forward` and `polar_backward` define reorientations,
while `azimuthal_forward` and `azimuthal_forward` are ignored.
"""
abstract type AbstractMotilityTwoStep <: AbstractMotility end

# just a wrapper to allow state to be mutable
mutable struct MotileState
    state
end
MotileState() = MotileState(TwoState())
@enum TwoState::Bool Forward Backward
Base.show(io::IO, ::MIME"text/plain", x::TwoState) = 
    x == Forward ? print(io, "Forward") : print(io, "Backward")
# choose at random between Forward and Backward if not specified
TwoState() = TwoState(Random.default_rng())
TwoState(rng::AbstractRNG) = TwoState(rand(rng, (true, false)))
# overload getproperty and setproperty! for more convenient access to state
function Base.getproperty(obj::AbstractMotilityTwoStep, sym::Symbol)
    if sym === :state
        return obj.motile_state.state
    else
        return getfield(obj, sym)
    end
end
function Base.setproperty!(value::AbstractMotilityTwoStep, name::Symbol, x)
    if name === :state
        return setfield!(value.motile_state, :state, x)
    else
        return setfield!(obj, name, x)
    end
end
# define rules for switching motile state
switch!(::AbstractMotilityOneStep) = nothing
switch!(m::AbstractMotilityTwoStep) = switch!(m.state)
switch!(s::TwoState) = (s = ~s; nothing)
Base.:~(x::TwoState) = TwoState(~Bool(x))

Base.@kwdef struct RunTumble <: AbstractMotilityOneStep
    speed = Degenerate(30.0)
    polar = Uniform(-π, π)
    azimuthal = Arccos(-1, 1)
end # struct

Base.@kwdef struct RunReverse <: AbstractMotilityTwoStep
    speed_forward = Degenerate(30.0)
    polar_forward = Degenerate(π)
    azimuthal_forward = Arccos(-1,1)
    speed_backward = speed_forward
    polar_backward = polar_forward
    azimuthal_backward = azimuthal_forward
    motile_state = MotileState()
end # struct

Base.@kwdef struct RunReverseFlick <: AbstractMotilityTwoStep
    speed_forward = Degenerate(30.0)
    polar_forward = Degenerate(π)
    azimuthal_forward = Arccos(-1,1)
    speed_backward = speed_forward
    polar_backward = [-π/2, π/2]
    azimuthal_backward = Arccos(-1,1)
    motile_state = MotileState()
end # struct 