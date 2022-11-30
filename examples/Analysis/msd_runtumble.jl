using BacteriaBasedModels
using Distributions
using Plots

θs = [π/6, π/4, π/3, π/2, π]

U = 30.0 # μm/s 
τ_run = 1.0 # s 
ω = 1 / τ_run

nmicrobes = 100
microbes = [
    [
        Microbe{3}(
            id = n, turn_rate = ω,
            motility = RunTumble(
                speed=Degenerate(U), polar=[θ,-θ],
                )
        ) for n in 1:nmicrobes 
    ] for θ in θs
]

dt = 0.05 # s 
L = 500.0 # μm
models = [
    initialise_model(;
        microbes = microbes[i],
        timestep = dt, extent = L
    ) for i in eachindex(microbes)
]

nsteps = round(Int, 100τ_run / dt)
adata = [:pos]
adfs = [run!(model, microbe_step!, nsteps; adata)[1] for model in models]

MSD = msd.(adfs; L=L)
ts = (1:nsteps).*dt

logslice = [1,2,5,10,25,50,100,250,500,1000]
plot(
    xlab = "Δt (s)",
    ylab = "MSD (μm²/s)",
    legend = :bottomright, legendtitle = "1-cosθ",
    scale = :log10
)
scatter!(ts[logslice], hcat(MSD...)[logslice,:],
    m=:x, ms=6, msw=2, lab=false)
for i in eachindex(θs)
    α = 1 - cos(θs[i])
    τ = τ_run / α
    D = U^2*τ / 3
    dr² = @. 2*U^2*τ^2 * (ts/τ - 1 + exp(-ts/τ))
    plot!(ts, dr², lab=round(α,digits=2), lc=i, lw=2)
end # for
plot!()