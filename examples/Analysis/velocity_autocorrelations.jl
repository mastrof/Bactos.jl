using Bactos
using Plots

U = 30.0 # μm/s
τ_run = 1.0 # s
ω = 1 / τ_run # 1/s
Δt = 0.01 # s
extent = 1e4 # μm

n = 200
microbes_runtumble = [
    Microbe{3}(id=i,
        turn_rate=ω, motility=RunTumble(speed=[U])
    )
    for i in 1:n
]
microbes_runrev = [
    Microbe{3}(id=n+i,
        turn_rate=ω, motility=RunReverse(speed_forward=[U])
    )
    for i in 1:n
]
microbes_runrevflick = [
    Microbe{3}(id=2n+i,
        turn_rate=ω, motility=RunReverseFlick(speed_forward=[U])
    )
    for i in 1:n
]

microbes = vcat(
    microbes_runtumble, microbes_runrev, microbes_runrevflick
)

model = initialise_model(;
    microbes,
    extent, periodic = true,
    timestep = Δt

)

nsteps = round(Int, 100τ_run / Δt)
adata = [:vel]
adf, = run!(model, nsteps; adata)

adf_runtumble = filter(:id => id -> model.agents[id].motility isa RunTumble, adf; view=true)
adf_runrev = filter(:id => id -> model.agents[id].motility isa RunReverse, adf; view=true)
adf_runrevflick = filter(:id => id -> model.agents[id].motility isa RunReverseFlick, adf; view=true)
adfs = [adf_runtumble, adf_runrev, adf_runrevflick]

Φ = hcat([autocorrelation(a,:vel) for a in adfs]...)

t = range(0, (nsteps-1)*Δt; step=Δt)
Φ_theoretical = hcat([
    exp.(-t ./ τ_run),
    exp.(-t ./ (τ_run / 2)),
    (1 .- t ./ (2τ_run)) .* exp.(-t ./ τ_run),
]...) # Taktikos et al. 2013 PLoS ONE

plot(
    xlims=(0,6τ_run), ylims=(-0.1, 1.05),
    xlab="Δt / τ",
    ylab="velocity autocorrelation",
)
plot!(t, Φ_theoretical, lw=2, lc=[1 2 3], label=["Run-Tumble" "Run-Reverse" "Run-Reverse-Flick"])
scatter!(t[1:15:end], Φ[1:15:end,:] ./ U^2, m=:x, mc=[1 2 3], label=false)
hline!([0.0], lw=0.8, ls=:dash, lc=:black, lab=false)