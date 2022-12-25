using Bactos
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

motility_fw = RunReverseFlick(motile_state = TwoStates(ForwardState()))
motility_bw = RunReverseFlick(motile_state = TwoStates(BackwardState()))
microbes = [
    XieNoisy{3}(id=1, turn_rate_forward=0, motility=motility_fw),
    XieNoisy{3}(id=2, turn_rate_backward=0, motility=motility_bw)
]

timestep = 0.1 # s
extent = 500.0 # μm
C₀ = 0.01 # μM
C₁ = 5.0-C₀ # μM 
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

_β(a) = motilestate(a.motility) == ForwardState() ? a.gain_forward : a.gain_backward
state(a::AbstractXie) = max(1 + _β(a)*a.state, 0)
adata = [state, :state_m, :state_z]
adf, = run!(model, microbe_step!, my_model_step!, nsteps; adata)
S = vectorize_adf_measurement(adf, :state)'
m = (vectorize_adf_measurement(adf, :state_m)')[:,1] # take only fw
z = (vectorize_adf_measurement(adf, :state_z)')[:,1] # take only fw

# response vs time for fw and bw modes
begin
    _green = palette(:default)[3]
    plot()
    x = (0:timestep:T) .- t₁
    plot!(
        x, S,
        lw=1.5, lab=["Forward" "Backward"]
    )
    plot!(ylims=(-0.1,4.5), ylab="Response", xlab="time (s)")
    plot!(twinx(),
        x, t -> concentration_field(t.+t₁,C₀,C₁,t₁,t₂),
        ls=:dash, lw=1.5, lc=_green, lab=false,
        tickfontcolor=_green,
        ylab="C (μM)", guidefontcolor=_green
    )
end

# methylation and dephosphorylation
begin
    x = (0:timestep:T) .- t₁
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