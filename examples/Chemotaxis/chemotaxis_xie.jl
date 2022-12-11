using BacteriaBasedModels
using Plots

θ(a,b) = a>b ? 1.0 : 0.0
function concentration_field(pos, model)
    C₀ = model.C₀
    C₁ = model.C₁
    t₁ = model.t₁
    t₂ = model.t₂
    dt = model.timestep
    t = model.t * dt
    concentration_field(t, C₀, C₁, t₁, t₂)
end
concentration_field(t,C₀,C₁,t₁,t₂) = C₀+C₁*θ(t,t₁)*(1-θ(t,t₂))
concentration_gradient(pos, model) = zero.(pos) # not used by Xie
concentration_time_derivative(pos, model) = 0.0 # not used by Xie

microbes = [XieNoisy{3}(id=1)]

timestep = 0.1 # s
extent = 500.0 # μm
C₀ = 0.01 # μM
C₁ = 2.5-C₀ # μM 
T = 60.0 # s
t₁ = 20.0 # s
t₂ = 40.0 # s
nsteps = round(Int, T/timestep)

model_properties = Dict(
    :compound_diffusivity => 608.0, # μm²/s
    :concentration_field => concentration_field,
    :concentration_gradient => concentration_gradient,
    :concentration_time_derivative => concentration_time_derivative,
    :C₀ => C₀,
    :C₁ => C₁,
    :t₁ => t₁,
    :t₂ => t₂,
    :t => 0,
)

model = initialise_model(;
    microbes, timestep,
    extent, model_properties
)

update_model!(model) = (model.t += 1)
my_model_step!(model) = model_step!(model; update_model!)

state(a::AbstractXie) = 1 + a.gain_forward*a.state
adata = [state, :state_m, :state_z]
adf, = run!(model, microbe_step!, my_model_step!, nsteps; adata)


# instantaneous tumbling vs time
begin
    _blue = palette(:default)[1]
    _orange = palette(:default)[2]
    ν₀ = microbes[1].turn_rate
    plot()
    x = (0:timestep:T) .- t₁
    y = max.(adf.state .* (ν₀), 0)
    plot!(
        x, y,
        lw=2, lab=false
    )
    plot!(ylims=(-0.1,4.5), ylab="ω (1/s)", xlab="time (s)")
    plot!(yguidefontcolor=_blue, ytickfontcolor=_blue)
    plot!(twinx(),
        x, t -> concentration_field(t.+t₁,C₀,C₁,t₁,t₂),
        ls=:dash, lw=1.5, lc=2, lab=false,
        tickfontcolor=_orange,
        ylab="C (μM)", guidefontcolor=_orange
    )
end

# methylation and dephosphorylation
begin
    x = (0:timestep:T) .- t₁
    m = adf.state_m
    z = adf.state_z
    τ_m = microbes[1].adaptation_time_m
    τ_z = microbes[1].adaptation_time_z
    M = m ./ τ_m
    Z = z ./ τ_z
    R = M .- Z
    plot(
        x, [M Z R],
        lw=2,
        lab=["m/τ_m" "z/τ_z" "m/τ_m - z/τ_z"],
        xlab="time (s)"
    )
end