using BacteriaBasedModels
using Plots

θ(a,b) = a>b ? 1. : 0.
conc_field(x,C,C₀,x₀) = C₀ + C*(1-θ(x,2x₀/3))*θ(x,x₀/3)

U = 30.0
motility = RunTumble(speed=Degenerate(U))

nmicrobes = 1
microbes = [
    CelaniNoisy{1}(id = i, gain=1.0, memory=0.5,
                   pos=(0.,),
                   vel=(U,), turn_rate=0.0,
                   motility=motility,
                   rotational_diffusivity=0.,
                   chemotactic_precision=0.0)
    for i in 1:nmicrobes
]


dt = 0.1 # s
L = 450.0 # μm
# field properties
C = 1.0 # μM 
C₀ = 0.0 # μm 
x₀ = y₀ = L # μm 
concentration_field(pos) = conc_field(pos[1], C, C₀, x₀)#*(1 + (rand()-0.5)/(5C))

model_properties = Dict(
    :concentration_field => (pos,_) -> concentration_field(pos),
    :concentration_gradient => (pos,_) -> 0.0,
    :concentration_time_derivative => (_,_) -> 0.0,
    :compound_diffusivity => 500.0
)

model = initialise_model(;
    microbes = microbes,
    timestep = dt,
    extent = L, periodic = true,
    model_properties = model_properties,
    random_positions = false,
)

state(a) = copy(a.state[4])
adata = [:pos, :vel, state]
nsteps = round(Int, L/(U*dt))-1
adf, = run!(model, microbe_step!, nsteps; adata)

traj = vectorize_adf_measurement(adf, :pos)
x = first.(traj)'
y = last.(traj)'

S = vectorize_adf_measurement(adf, :state)'
c = concentration_field.([(x,) for x in range(0,L;length=nsteps+1)])./C
plot((0:nsteps).*dt, c)
plot!((0:nsteps).*dt, 1 .- S, m=:c, msw=0)