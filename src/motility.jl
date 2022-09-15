export
    AbstractMotility, AbstractMotilityOneStep, AbstractMotilityTwoStep,
    RunTumble, RunReverse, RunReverseFlick

abstract type AbstractMotility end
abstract type AbstractMotilityOneStep end
abstract type AbstractMotilityTwoStep end

Base.@kwdef struct RunTumble <: AbstractMotilityOneStep
    speed = Degenerate(1.0)
    angle = Uniform(-π, π)
end # struct

Base.@kwdef struct RunReverse <: AbstractMotilityOneStep
    speed = Degenerate(1.0)
    angle = Degenerate(π)
end # strut 

Base.@kwdef struct RunReverseFlick <: AbstractMotilityTwoStep
    speed = Degenerate(1.0)
    angle_reverse = Degenerate(π)
    angle_flick = Degenerate(π/2)
    motile_state::Vector{Int} = rand(0:1, 1)
end # struct 