module Bactos

using Agents, Agents.Pathfinding
using CellListMap.PeriodicSystems
using Distributions
using LinearAlgebra
using StaticArrays
using Random
using Rotations

using OrdinaryDiffEq: ODEProblem, DEIntegrator, init, Tsit5, step!
using FiniteDifferences: central_fdm

# Core structures
include("distributions.jl")
include("motility.jl")
include("microbes.jl")

# Utility routines
include("utils.jl")

# ABM setup
include("model.jl")

# Extras
include("diffeq.jl")
include("pathfinder.jl")
include("obstacles_spheres.jl")
include("celllistmap.jl")
include("encounters.jl")

# Stepping
include("rotations.jl")
include("step_microbes.jl")
include("finite_differences.jl")

# Chemotaxis models
include("brown-berg.jl")
include("brumley.jl")
include("celani.jl")
include("xie.jl")

# Measurements
include("drift.jl")
include("msd.jl")
include("correlation_functions.jl")
end # module