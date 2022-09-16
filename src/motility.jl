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
One-step motility patterns (such as run-tumble and run-reverse).
Subtypes must have at least the following fields:
- `speed`: distribution of microbe speed, new values extracted after each turn
- `yaw`: distribution of in-plane reorientations
- `pitch`: distribution of out-of-plane reorientations
For 2-dimensional microbe types, only `yaw` defines reorientations.
"""
abstract type AbstractMotilityOneStep end
abstract type AbstractMotilityTwoStep end

Base.@kwdef struct RunTumble <: AbstractMotilityOneStep
    speed = Degenerate(1.0)
    yaw = Uniform(-π, π) # in-plane
    pitch = Uniform(0,2π) # out-of-plane
end # struct

Base.@kwdef struct RunReverse <: AbstractMotilityOneStep
    speed = Degenerate(1.0)
    yaw = Degenerate(π)
    pitch = Uniform(0,2π)
end # strut 

"""
    AbstractMotilityTwoStep
Two-step motility patterns (such as run-reverse-flick).
Subtypes must have at least the following fields:
- `speed`: distribution of microbe speed, new values extracted after each turn
- `yaw_0`: distribution of in-plane reorientations for motile state 0
- `pitch_0`: distribution of out-of-plane reorientations for motile state 0
- `yaw_1`: distribution of in-plane reorientations for motile state 1
- `pitch_1`: distribution of out-of-plane reorientations for motile state 1
- `motile_state::Vector{Int}`: defines current motile state (`[0]` or `[1]`)
For 2-dimensional microbe types, only `yaw_0` and `yaw_1` define reorientations.
"""
Base.@kwdef struct RunReverseFlick <: AbstractMotilityTwoStep
    speed = Degenerate(1.0)
    yaw_0 = Degenerate(π)
    pitch_0 = Uniform(0,2π)
    yaw_1 = Degenerate(π/2)
    pitch_1 = Uniform(0,2π)
    motile_state::Vector{Int} = rand(0:1, 1)
end # struct 