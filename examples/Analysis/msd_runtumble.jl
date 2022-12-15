using BacteriaBasedModels
using Distributions
using Plots

θs = [π/6, π/4, π/3, π/2, π]
α = cos.(θs)

U = 30.0 # μm/s 
τ = 1.0 # s 
ω = 1 / τ

nmicrobes = 100
microbes = [
    [
        Microbe{3}(
            id = n, turn_rate = ω,
            motility = RunTumble(speed=[U], polar=[θ,-θ])
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

MSD = hcat(msd.(adfs; L=L)...)

begin  
    t = (1:nsteps).*dt
    T = @. τ / (1-α')
    s = t ./ T
    D = @. U^2*T/3
    MSD_theoretical = @. 6D*T * (s - 1 + exp(-s))
    logslice = round.(Int, exp10.(range(0,3,length=10)))
    plot(
        xlab = "Δt / τ",
        ylab = "MSD / (Uτ)²",
        legend = :bottomright, legendtitle = "1-α",
        scale = :log10,
        yticks = exp10.(-2:2:2),
        xticks = exp10.(-2:2)
    )
    scatter!(t[logslice,:]./τ, MSD[logslice,:]./(U*τ)^2,
        m=:x, ms=6, msw=2, lab=false, lc=axes(α,1)'
    )
    plot!(t./τ, MSD_theoretical./(U*τ)^2,
        lw=2, lab=round.(1 .- α,digits=2)', lc=axes(α,1)'
    )
end