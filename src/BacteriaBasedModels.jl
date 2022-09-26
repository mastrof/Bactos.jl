module BacteriaBasedModels

using Agents
using Distributions
using LinearAlgebra
using StaticArrays
using Random
using Rotations

using OrdinaryDiffEq: ODEProblem, DEIntegrator, init, Tsit5, step!
using FiniteDifferences: central_fdm

# Utility routines
include("utils.jl")

# Core structures
include("distributions.jl")
include("motility.jl")
include("microbes.jl")

# ABM setup
include("model.jl")

# Stepping
export run! # from Agents
include("rotations.jl")
include("step_microbes.jl")
include("step_model.jl")
include("finite_differences.jl")

# Chemotaxis models
include("brown-berg.jl")
include("brumley.jl")

# Measurements
include("drift.jl")
include("msd.jl")
include("correlation_functions.jl")
end # module