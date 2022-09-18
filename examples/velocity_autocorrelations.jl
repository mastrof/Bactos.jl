using BacteriaBasedModels
using Plots

U = 30.0 # μm/s
τ_run = 1.0 # s
ω = 1 / τ_run # 1/s
Δt = 0.1 # s
L = 1e4 # μm

n = 100
microbes_runtumble = [
    Microbe{3}(id=i,
        turn_rate=ω, vel=rand_vel(3).*U,
        motility=RunTumble(speed=Degenerate(U))
    )
    for i in 1:n
]
microbes_runrev = [
    Microbe{3}(id=n+i,
        turn_rate=ω, vel=rand_vel(3).*U,
        motility=RunReverse(speed=Degenerate(U))
    )
    for i in 1:n
]
microbes_runrevflick = [
    Microbe{3}(id=2n+i,
        turn_rate=ω, vel=rand_vel(3).*U,
        motility=RunReverseFlick(speed=Degenerate(U))
    )
    for i in 1:n
]

microbes = vcat(
    microbes_runtumble, microbes_runrev, microbes_runrevflick
)

model = initialise_model(;
    microbes = microbes,
    timestep = Δt,
    extent = L, periodic = true
)

nsteps = 10_000
adata = [:pos, :vel]
adf, = run!(model, microbe_step!, nsteps; adata)

adf_runtumble = filter(:id => id -> 1≤id≤n, adf; view=true)
adf_runrev = filter(:id => id -> n+1≤id≤2n, adf; view=true)
adf_runrevflick = filter(:id => id -> 2n+1≤id≤3n, adf; view=true)

adfs = [adf_runtumble, adf_runrev, adf_runrevflick]

Φ = hcat([autocorrelation(a,:vel) for a in adfs]...)

t = range(0, (nsteps-1)*Δt; step=Δt)
ϕ = hcat([
    exp.(-t ./ τ_run),
    exp.(-t ./ (τ_run / 2)),
    (1 .- t ./ (2τ_run)) .* exp.(-t ./ τ_run),
]...) # Taktikos et al. 2013 PLoS ONE

plot(
    xlims=(0,6τ_run), ylims=(-0.1, 1.05), legend=false,
    xlab="Δt / τ",
    ylab="velocity autocorrelation",
)
plot!(t, ϕ, lw=2, lc=[1 2 3])
scatter!(t[1:2:end], Φ[1:2:end,:] ./ U^2, m=:x, mc=[1 2 3])