using Bactos
using Plots

L = 50.0
τ_run = 1.0
ω = 1 / τ_run
nmicrobes = 6
microbes = [
    Microbe{2}(
        id=i, vel=rand_vel(2), turn_rate=ω,
        motility=RunTumble()
    ) for i in 1:nmicrobes
]

model = initialise_model(;
    microbes = microbes,
    timestep = 1.0,
    extent = L, periodic = false,
    random_positions = true,
)

nsteps = 200
adata = [:pos]
adf, = run!(model, microbe_step!, nsteps; adata)

trajectories = vectorize_adf_measurement(adf, :pos)
x = first.(trajectories)
y = last.(trajectories)

plot(x', y',
    lw=0.5,
    lab=false, ticks=false, lims=(0,L)
)