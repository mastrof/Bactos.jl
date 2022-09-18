using BacteriaBasedModels
using Agents: run!
using Plots

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


L = 50.0
vel = (1.0,)
τ_run = 1.0
ω = 1 / τ_run
m = Microbe{1}(id=1, vel=vel, turn_rate=ω)

Δt = 0.1
model = initialise_model(;
    microbes = [m],
    timestep = Δt,
    extent = L, periodic = false,
)

nsteps = 10_000
adata = [:pos, :vel]
adf, = run!(model, microbe_step!, nsteps; adata)

x = adf[!,:pos] .|> first
plot(x.-x[1], lab=false, xlab="steps", ylab="Δx")