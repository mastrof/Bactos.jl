export AbstractMicrobe, Microbe

"""
    YourMicrobeType{D} <: AbstractMicrobe{D} where {D<:Integer}
All microbe types in BacteriaBasedModels.jl simulations must be instances
of user-defined types that are subtypes of `AbstractMicrobe`.
The parameter `D` defines the dimensionality of the space in which the
microbe type lives (1, 2 and 3 are currently supported).

All microbe types should have the following fields:
`id::Int` → id of the microbe
`pos::NTuple{D,Float64}` → position of the microbe
`vel::NTuple{D,Float64}` → velocity of the microbe
`motility` → motile pattern of the microbe
`turn_rate::Float64` → average reorientation rate of the microbe
`rotational_diffusivity::Float64` → rotational diffusion coefficient
"""
abstract type AbstractMicrobe{D} <: AbstractAgent where {D<:Integer} end


Base.@kwdef mutable struct Microbe{D} <: AbstractMicrobe{D}
    id::Int
    pos::NTuple{D,Float64} = ntuple(_ -> 0.0, D)
    vel::NTuple{D,Float64} = rand_vel(D)
    motility = RunTumble()
    turn_rate::Float64 = 1.0
    state::Float64 = 0.0
    rotational_diffusivity::Float64 = 0.0
    radius::Float64 = 0.0
end # struct 