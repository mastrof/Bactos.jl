using BacteriaBasedModels
using Plots

conc_field(x,y,C,σ,x₀,y₀) = C * exp(-((x-x₀)^2+(y-y₀)^2)/(2σ^2))
function conc_grad(x,y,C,σ,x₀,y₀)
    cfield = conc_field(x,y,C,σ,x₀,y₀)
    σ² = σ*σ
    return [
        -(x-x₀)/σ² * cfield,
        -(y-y₀)/σ² * cfield
    ]
end # function

nmicrobes = 10
microbes = [MicrobeBrumley{2}(id=i) for i in 1:nmicrobes]

dt = 0.1 # s
L = 500.0 # μm
# field properties
C = 0.1 # μM 
σ = 20.0 # μm 
x₀ = y₀ = L/2 # μm 
concentration_field(pos) = conc_field(pos[1], pos[2], C, σ, x₀, y₀)
concentration_gradient(pos) = conc_grad(pos[1], pos[2], C, σ, x₀, y₀)

model_properties = Dict(
    :concentration_field => (pos,_) -> concentration_field(pos),
    :concentration_gradient => (pos,_) -> concentration_gradient(pos),
    :concentration_time_derivative => (_,_) -> 0.0,
    :compound_diffusivity => 500.0, # μm²/s
)

model = initialise_model(;
    microbes = microbes,
    timestep = dt,
    extent = L, periodic = false,
    model_properties = model_properties
)

my_microbe_step!(microbe, model) = microbe_step!(
    microbe, model;
    affect! = brumley_affect!,
    turnrate = brumley_turnrate
)

adata = [:pos]
nsteps = 1000
adf, = run!(model, my_microbe_step!, nsteps; adata)

traj = vectorize_adf_measurement(adf, :pos)
x = first.(traj)'
y = last.(traj)'
contourf(
    0:L/50:L, 0:L/50:L,
    (x,y) -> concentration_field((x,y)),
    color=:bone, ratio=1
)
plot!(
    x, y, lw=0.5,
    colorbar=false, legend=false,
    xlims=(0,L), ylims=(0,L)
)