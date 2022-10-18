using BacteriaBasedModels
using Plots

conc_field(x,y,C,σ,x₀,y₀) = C * exp(-((x-x₀)^2+(y-y₀)^2)/(2σ^2))
function conc_grad(x,y,C,σ,x₀,y₀)
    cfield = conc_field(x,y,C,σ,x₀,y₀)
    σ² = σ*σ
    return [
        0.0,
        0.0,#-(x-x₀)/σ² * cfield,
        0.0#-(y-y₀)/σ² * cfield
    ]
end # function

using Random
Random.seed!(3)

U = 30.0
motility = RunTumble(speed=Degenerate(U))

nmicrobes = 100
microbes = [
    Celani{2}(id = i, gain=50.0, memory=0.5,
                   vel=rand_vel(2).*U,
                   motility=motility,
                   rotational_diffusivity=0.1,)
                   #chemotactic_precision=0.0)
    for i in 1:nmicrobes
]


dt = 0.1 # s
L = 1500.0 # μm
# field properties
C = 50.0 # μM 
σ = 35.0 # μm 
x₀ = y₀ = L/2 # μm 
concentration_field(pos) = conc_field(pos[1], pos[2], C, σ, x₀, y₀)
concentration_gradient(pos) = conc_grad(pos[1], pos[2], C, σ, x₀, y₀)

model_properties = Dict(
    :concentration_field => (pos,_) -> concentration_field(pos),
    :concentration_gradient => (pos,_) -> concentration_gradient(pos),
    :concentration_time_derivative => (_,_) -> 0.0,
    :compound_diffusivity => 500.0
)

model = initialise_model(;
    microbes = microbes,
    timestep = dt,
    extent = L, periodic = false,
    model_properties = model_properties
)

state(a) = copy(a.state)
adata = [:pos, :vel, state]
nsteps = 1000
adf, = run!(model, microbe_step!, nsteps; adata)

traj = vectorize_adf_measurement(adf, :pos)
x = first.(traj)'
y = last.(traj)'


contourf(
    0:L/50:L, 0:L/50:L,
    (x,y) -> concentration_field((x,y)),
    color=:bone, ratio=1
)
plot!(
    x, y, lw=0.4,
    colorbar=false, legend=false,
    xlims=(0,L), ylims=(0,L),
    bgcolor=:black, axis=false
)
