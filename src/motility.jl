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
- `polar`: distribution of in-plane reorientations
- `azimuthal`: distribution of out-of-plane reorientations
For 2-dimensional microbe types, only `polar` defines reorientations.
"""
abstract type AbstractMotilityOneStep <: AbstractMotility end

"""
    AbstractMotilityTwoStep
Two-step motility patterns (such as run-reverse-flick).
Subtypes must have at least the following fields:
- `speed`: distribution of microbe speed, new values extracted after each turn
- `polar_0`: distribution of in-plane reorientations for motile state 0
- `azimuthal_0`: distribution of out-of-plane reorientations for motile state 0
- `polar_1`: distribution of in-plane reorientations for motile state 1
- `azimuthal_1`: distribution of out-of-plane reorientations for motile state 1
- `motile_state::Vector{Int}`: defines current motile state (`[0]` or `[1]`)
For 2-dimensional microbe types, only `polar_0` and `polar_1` define reorientations.
"""
abstract type AbstractMotilityTwoStep <: AbstractMotility end

Base.@kwdef struct RunTumble <: AbstractMotilityOneStep
    speed = Degenerate(1.0)
    polar = Uniform(-π, π) # in-plane
    azimuthal = Arccos(-1, 1) # out-of-plane
end # struct

Base.@kwdef struct RunReverse <: AbstractMotilityOneStep
    speed = Degenerate(1.0)
    polar = Degenerate(π)
    azimuthal = Arccos(-1, 1)
end # strut 

Base.@kwdef struct RunReverseFlick <: AbstractMotilityTwoStep
    speed = Degenerate(1.0)
    polar_0 = Degenerate(π)
    azimuthal_0 = Arccos(-1, 1)
    polar_1 = Degenerate(π/2)
    azimuthal_1 = Arccos(-1, 1)
    motile_state::Vector{Int} = rand(0:1, 1)
end # struct 