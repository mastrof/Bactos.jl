module BacteriaBasedModels

using Agents
using Distributions
using LinearAlgebra
using StaticArrays
using Random
using Rotations

# Utility routines
include("utils.jl")

# Core structures
include("distributions.jl")
include("motility.jl")
include("microbes.jl")

# ABM setup
include("model.jl")

# Stepping
include("step_microbes.jl")
end # module