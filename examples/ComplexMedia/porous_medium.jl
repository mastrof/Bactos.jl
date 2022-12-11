using BacteriaBasedModels
using DelimitedFiles
using StaticArrays
using LinearAlgebra
using Plots

# overlap between spheres
isoverlapping(p1, p2, r1, r2) = (norm(p1 .- p2) ≤ r1+r2)
isoverlapping(a, b) = isoverlapping(a.pos, b.pos, a.radius, b.radius)

# draw a circle
function circleShape(x₀,y₀,r,n=500)
    θ = LinRange(0, 2π, n)
    x₀ .+ r.*sin.(θ), y₀ .+ r.*cos.(θ)
end # function


# Physical parameters
timestep = 0.1 # s
extent = (1000.0, 500.0) # μm
periodic = false
microbe_radius = 0.5 # μm
ω = 1.0 # 1/s
U = 30.0 # μm/s 
motility = RunTumble(speed = Degenerate(U))
Drot = 0.1 # rad²/s
n_microbes = 6

# Initialise obstacles (read configuration from file) 
obstacle_data = readdlm("phi04_rmin10_Lx1000_Ly500.dat")
bodyrad = obstacle_data[:,1] # μm
max_radius = maximum(bodyrad)
bodypos = [Tuple(obstacle_data[i,2:3]) for i in axes(obstacle_data,1)] # μm
bodies = [
    ObstacleSphere(pos, r, glide!) for (r,pos) in zip(bodyrad,bodypos)
]

# Initialise microbes
microbes = [
    Microbe{2}(
        id=i, pos=Tuple(rand(2).*extent), vel=rand_vel(2).*U,
        turn_rate=1.0, radius=0.5,
        motility=RunTumble(speed=Degenerate(U)),
        rotational_diffusivity=Drot
    ) for i in 1:n_microbes
]
# Update microbe positions to avoid overlap with obstacles
for m in microbes
    while any(map(b -> isoverlapping(m,b), bodies))
        m.pos = Tuple(rand(2) .* extent)
    end # while
end # for

# Initialise neighbor list
cutoff_radius = 2 * (max_radius + microbe_radius + U*timestep)
neighborlist = init_neighborlist(microbes, bodies, extent, periodic, cutoff_radius)

model_properties = Dict(
    :bodies => bodies,
    :neighborlist => neighborlist
)

model = initialise_model(;
    microbes, timestep, extent, periodic, model_properties,
    random_positions = false
)

# required to update positions in the neighborlist
function affect!(microbe, model)
    update_microbepos_neighborlist!(microbe, model)
end # function

my_microbe_step!(microbe, model) = microbe_step!(
    microbe, model; affect!
)

function update_model!(model)
    surface_interaction!(model)
end # function

my_model_step!(model) = model_step!(model; update_model!)


adata = [:pos, :vel]
adf, = run!(model, my_microbe_step!, my_model_step!, 2000; adata)

traj = vectorize_adf_measurement(adf, :pos)
x = first.(traj)'
y = last.(traj)'

plot(
    xlims=(0,extent[1]), ylims=(0,extent[2]),
    palette=:Dark2, legend=false,
    bgcolor=:black, grid=false, axis=false,
    ratio=1)

for body in bodies
    plot!(
        circleShape(body.pos..., body.radius),
        seriestype=:shape, lw=0, lab=false,
        c=:white, fillalpha=0.25,
    )
end # for

plot!(x,y)