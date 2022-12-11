export
    AbstractMotility, AbstractMotilityOneStep, AbstractMotilityTwoStep,
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
- `motile_state`: defines current motile state (`ForwardState` or `BackwardState`)
For 2-dimensional microbe types, only `polar_forward` and `polar_backward` define reorientations,
while `azimuthal_forward` and `azimuthal_forward` are ignored.
"""
abstract type AbstractMotilityTwoStep <: AbstractMotility end

abstract type AbstractMotileState end
struct ForwardState; end
struct BackwardState; end
ForwardBackward = Union{ForwardState,BackwardState}
Base.@kwdef mutable struct TwoStates <: AbstractMotileState
    state::ForwardBackward = rand((ForwardState, BackwardState))()
end
switch(::ForwardState) = BackwardState()
switch(::BackwardState) = ForwardState()
switch!(s::TwoStates) = (s.state = switch(s.state))

Base.@kwdef struct RunTumble <: AbstractMotilityOneStep
    speed = Degenerate(30.0)
    polar = Uniform(-π, π) # in-plane
    azimuthal = Arccos(-1, 1) # out-of-plane
end # struct

#=
Base.@kwdef struct RunReverse <: AbstractMotilityOneStep
    speed = Degenerate(1.0)
    polar = Degenerate(π)
    azimuthal = Arccos(-1, 1)
end # strut 
=#

Base.@kwdef struct RunReverse <: AbstractMotilityTwoStep
    speed_forward = Degenerate(30.0)
    polar_forward = Degenerate(π)
    azimuthal_forward = Arccos(-1,1)
    speed_backward = speed_forward
    polar_backward = polar_forward
    azimuthal_backward = azimuthal_forward
    motile_state = TwoStates()
end # struct

Base.@kwdef struct RunReverseFlick <: AbstractMotilityTwoStep
    speed_forward = Degenerate(30.0)
    polar_forward = Degenerate(π)
    azimuthal_forward = Arccos(-1,1)
    speed_backward = speed_forward
    polar_backward = Degenerate(π/2)
    azimuthal_backward = Arccos(-1,1)
    motile_state = TwoStates()
end # struct 