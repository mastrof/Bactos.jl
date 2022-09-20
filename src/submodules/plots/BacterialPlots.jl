export BacterialPlots

module BacterialPlots

using BacteriaBasedModels
using Plots
using LinearAlgebra: norm

# Default plot properties
default(
    thickness_scaling = 1.5,
    guidefontsize = 12,
    tickfontsize = 12,
    legendfontsize = 8,
    grid = false,
    framestyle = :box,
    minorticks = true,
    tick_direction = :in,
    color_palette = :Dark2,
    margin = 3.0Plots.mm
)

include("bacteria_2d.jl")

end # module