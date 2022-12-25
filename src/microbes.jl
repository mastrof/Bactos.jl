export AbstractMicrobe, Microbe

"""
    AbstractMicrobe{D} <: AbstractAgent where {D<:Integer}
    YourMicrobeType{D} <: AbstractMicrobe{D}
All microbe types in Bactos.jl simulations must be instances
of user-defined types that are subtypes of `AbstractMicrobe`.
The parameter `D` defines the dimensionality of the space in which the
microbe type lives (1, 2 and 3 are currently supported).

All microbe types should have the following fields:
`- id::Int` → id of the microbe
`- pos::NTuple{D,Float64}` → position of the microbe
`- vel::NTuple{D,Float64}` → velocity of the microbe
`- motility` → motile pattern of the microbe
`- turn_rate::Float64` → average reorientation rate of the microbe
`- rotational_diffusivity::Float64` → rotational diffusion coefficient
"""
abstract type AbstractMicrobe{D} <: AbstractAgent where {D<:Integer} end


"""
    Microbe{D} <: AbstractMicrobe{D}
Basic microbe type for simple simulations.

Default parameters:
- `id::Int` → identifier used internally by Agents.jl
- `pos::NTuple{D,Float64} = ntuple(zero,D)` → position
- `motility = RunTumble()` → motile pattern
- `vel::NTuple{D,Float64} = rand_vel(D) .* rand(motility.speed)` → velocity vector
- `turn_rate::Float64 = 1.0` → frequency of reorientations
- `state::Float64` → generic variable for a scalar internal state
- `rotational_diffusivity::Float64 = 0.0` → rotational diffusion coefficient
- `radius::Float64 = 0.0` → equivalent spherical radius of the microbe
"""
Base.@kwdef mutable struct Microbe{D} <: AbstractMicrobe{D}
    id::Int
    pos::NTuple{D,Float64} = ntuple(zero, D)
    motility = RunTumble()
    vel::NTuple{D,Float64} = rand_vel(D, motility)
    turn_rate::Float64 = 1.0
    state::Float64 = 0.0
    rotational_diffusivity::Float64 = 0.0
    radius::Float64 = 0.0
end # struct 