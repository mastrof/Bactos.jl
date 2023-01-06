using Bactos
using DelimitedFiles
using LinearAlgebra
using Plots

# overlap between spheres
isoverlapping(p1, p2, r1, r2) = (norm(p1 .- p2) ≤ r1+r2)
isoverlapping(a, b) = isoverlapping(a.pos, b.pos, a.radius, b.radius)

# draw a circle
function circleShape(x₀,y₀,r,n=50)
    θ = LinRange(0, 2π, n)
    x₀ .+ r.*sin.(θ), y₀ .+ r.*cos.(θ)
end # function


# Physical parameters
timestep = 0.1 # s
extent = (1000.0, 500.0) # μm
periodic = false
n_microbes = 6

# Initialise obstacles (read configuration from file) 
obstacle_data = readdlm("phi065_rmin5_Lx1000_Ly500.dat")
bodyrad = obstacle_data[:,1] # μm
max_radius = maximum(bodyrad)
bodypos = [Tuple(obstacle_data[i,2:3]) for i in axes(obstacle_data,1)] # μm
bodies = [
    ObstacleSphere(pos, r) for (r,pos) in zip(bodyrad,bodypos)
]

# Initialise microbes at x=0
microbes = [Celani{2}(
    id=i, pos=(0,rand()*extent[2])) for i in 1:n_microbes
]
# Update microbe positions to avoid overlap with obstacles
for m in microbes
    while any(map(b -> isoverlapping(m,b), bodies))
        m.pos = (0, rand()*extent[2])
    end # while
end # for

# Setup concentration field
C₀=0.0
C₁=10.0
concentration_field(x,y,C₀,C₁,Lx) = C₀+(C₁-C₀)/Lx*x
function concentration_field(pos, model)
    C₀, C₁ = model.cfield_params
    Lx = model.space.extent[1]
    x, y = pos
    concentration_field(x,y,C₀,C₁,Lx)
end
function concentration_gradient(pos,model)
    C₀, C₁ = model.cfield_params
    Lx = model.space.extent[1]
    return ((C₁-C₀)/Lx, 0.0)
end

model_properties = Dict(
    :cfield_params => (C₀, C₁),
    :concentration_field => concentration_field,
    :concentration_gradient => concentration_gradient,
)

model = initialise_model(;
    microbes, timestep, extent, periodic, model_properties,
    random_positions = false
)
add_pathfinder!(model, 0.5, bodies)

adata = [:pos]
@time adf, = run!(model, 8000; adata)

traj = vectorize_adf_measurement(adf, :pos)
x = first.(traj)'
y = last.(traj)'

plot(
    xlims=(0,extent[1]), ylims=(0,extent[2]),
    palette=:Dark2, legend=false,
    bgcolor=:black, grid=false, axis=false,
    colorbar=:bottom, colorbartitle="C (μM)",
    ratio=1)

contourf!(
    0:extent[1], 0:extent[2],
    (x,y) -> concentration_field(x,y,C₀,C₁,extent[1]),
    color=:cividis, levels=100
)

for body in bodies
    plot!(
        circleShape(body.pos..., body.radius),
        seriestype=:shape, lw=0, lab=false,
        c=:black, fillalpha=0.5,
    )
end # for

plot!(x,y, lc=(1:n_microbes)')
scatter!(x[end:end,:], y[end:end,:], mc=(1:n_microbes)', ms=8, msc=:black)