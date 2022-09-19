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
L = 100.0
nmicrobes = 8
microbes = [Microbe{1}(id=i, pos=(L/2,)) for i in 1:nmicrobes]

model = initialise_model(;
    microbes = microbes,
    timestep = dt,
    extent = L, periodic = false,
    random_positions = false
)

nsteps = 1000
adata = [:pos]
adf, = run!(model, microbe_step!, nsteps; adata)

x = vectorize_adf_measurement(adf, :pos) .|> first
plot(
    (0:nsteps).*dt, x',
    legend = false,
    xlab = "time",
    ylab = "position"
)