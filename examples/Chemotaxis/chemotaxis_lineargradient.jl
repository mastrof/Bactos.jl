using BacteriaBasedModels
using Plots

#== CONCENTRATION FIELD ==#
conc_field(x,y,C₀,∇C) = C₀ + ∇C*x
function conc_grad(x,y,C₀,∇C)
    return [∇C, 0.0]
end # function

#== DOMAIN PARAMETERS ==#
Lx, Ly = 1000.0, 500.0 # μm
extent = (Lx, Ly) # μm
periodic = false

#== MICROBE POPULATION ==#
n = 50
microbes_brumley = [
    MicrobeBrumley{2}(id=i, pos=(0,rand()*Ly), chemotactic_precision=1)
    for i in 1:n
]
microbes_brown = [
    MicrobeBrownBerg{2}(id=n+i, pos=(0,rand()*Ly))
    for i in 1:n
]
MicrobeType = Union{MicrobeBrumley{2}, MicrobeBrownBerg{2}}
microbes = MicrobeType[microbes_brumley; microbes_brown]

#== FIELD PARAMETERS ==#
C₀ = 0.0 # μM 
∇C = 0.01 # μM/μm
concentration_field(pos) = conc_field(pos..., C₀, ∇C)
concentration_gradient(pos) = conc_grad(pos..., C₀, ∇C)

#== MODEL SETUP ==#
timestep = 0.1 # s
model_properties = Dict(
    :concentration_field => (pos,_) -> concentration_field(pos),
    :concentration_gradient => (pos,_) -> concentration_gradient(pos),
    :concentration_time_derivative => (_,_) -> 0.0,
    :compound_diffusivity => 500.0, # μm²/s
)

model = initialise_model(;
    microbes,
    timestep,
    extent, periodic,
    model_properties,
    random_positions = false
)

#== RUN! ==#
adata = [:pos]
nsteps = 1000
adf, = run!(model, microbe_step!, nsteps; adata)

#== POST PROCESSING ==#
traj = vectorize_adf_measurement(adf, :pos)
x = first.(traj)'
y = last.(traj)'

#== PLOTTING ==#
lc = [typeof(m) <: MicrobeBrumley ? 1 : 2 for m in microbes] |> permutedims

for t in axes(x,1)[1:4:end]
    contourf(
        0:Lx/50:Lx, 0:Ly/50:Ly,
        (x,y) -> concentration_field((x,y)),
        color=:bone, ratio=1,
        xlims=(0,Lx), ylims=(0,Ly),
        xlab="x (μm)", ylab="y (μm)",
        colorbar_title="C (μM)", levels=100,
    )

    t₀ = max(1, t-20)
    xt = @view x[t₀:t,:]
    yt = @view y[t₀:t,:]
    plot!(
        xt, yt,
        lw=0.5, lc=lc,
        legend=false,
    )

    it = lpad(t, 4, '0')
    savefig("frame_$it.png")
end