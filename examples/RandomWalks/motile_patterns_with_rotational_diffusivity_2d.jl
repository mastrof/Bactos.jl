using Bactos
using Plots

L = 800.0
τ_run = 1.0
ω = 1 / τ_run
U = 30.0
D_rot = 0.02

n = 1

m_rt = [
    Microbe{2}(
        id=i, turn_rate=ω,
        motility=RunTumble(speed=[U]),
        rotational_diffusivity=D_rot
    ) for i in 1:n
]

m_rr = [
    Microbe{2}(
        id=n+i, turn_rate=ω,
        motility=RunReverse(speed_forward=[U]),
        rotational_diffusivity=D_rot
    ) for i in 1:n
]

m_rrf = [
    Microbe{2}(
        id=2n+i, turn_rate=ω,
        motility=RunReverseFlick(speed_forward=[U]),
        rotational_diffusivity=D_rot
    ) for i in 1:n
]

microbes = vcat(m_rt, m_rr, m_rrf)

Δt = τ_run / 10
model = initialise_model(;
    microbes = microbes,
    timestep = Δt,
    extent = L, periodic = false,
)

nsteps = 500
adata = [:pos]
adf, = run!(model, nsteps; adata)

trajectories = vectorize_adf_measurement(adf, :pos)
x = first.(trajectories)
y = last.(trajectories)

plot(x', y',
    lw=0.5,
    ticks=false, lims=(0,L),
    lab=["Run-Tumble" "Run-Reverse" "Run-Reverse-Flick"],
    legend_foreground_color=:white, legendfontsize=5
)