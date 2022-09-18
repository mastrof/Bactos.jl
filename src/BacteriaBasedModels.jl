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
include("rotations.jl")
include("step_microbes.jl")

# Measurements
include("msd.jl")
include("correlation_functions.jl")
end # module