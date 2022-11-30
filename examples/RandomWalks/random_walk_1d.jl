using BacteriaBasedModels
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

dt = 0.1
L = 1e6
nmicrobes = 8
microbes = [Microbe{1}(id=i, pos=(L/2,)) for i in 1:nmicrobes]

model = initialise_model(;
    microbes = microbes,
    timestep = dt,
    extent = L, periodic = false,
)

nsteps = 1000
adata = [:pos]
adf, = run!(model, microbe_step!, nsteps; adata)

x = first.(vectorize_adf_measurement(adf, :pos))'
x₀ = x[1:1,:]
Δx = x .- x₀
plot(
    (0:nsteps).*dt, Δx,
    legend = false,
    xlab = "time",
    ylab = "position"
)