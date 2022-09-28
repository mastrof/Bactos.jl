using BacteriaBasedModels
using Agents: get_spatial_index
using OrdinaryDiffEq: get_du!
using Plots
default(
    thickness_scaling = 1.5,
    guidefontsize = 12,
    tickfontsize = 12,
    legendfontsize = 8,
    grid = false,
    framestyle = :box,
    minorticks = true,
    tick_direction = :in,
    color_palette = :Dark2,
    margin = 3.0Plots.mm
)

nmicrobes = 20
microbes = [MicrobeBrownBerg{1}(id=i) for i in 1:nmicrobes]

extent = 1000.0
spacing = 0.5
xs = range(-2*spacing, extent+2*spacing; step=spacing)
timestep = 0.1

x₀ = extent/2
C = 10.0
σ = 10.0
D = 10.0
β = 0.004
u₀ = @. C * exp(-(xs-x₀)^2 / (2*σ^2))
∇u = zero(u₀) # to be used as storage in model
∂ₜu = zero(u₀) # to be used as storage in model
finitediff!(∇u, u₀, 1/spacing)

function odestep!(du, u, p, t)
    β, D, _dx, = p
    a = D * _dx * _dx
    # diffusion
    laplacian!(du, u, a)
    # decay
    @. du -= β*u
    # absorbing walls
    du[1] = du[2] = du[end] = du[end-1] = 0.0
end # function

nsteps = round(Int, 500 / timestep)
ode_integrator = initialise_ode(odestep!, u₀, (β, D, 1/spacing);
                                dtmax = spacing^2/2D,
                                saveat = (0:nsteps) .* timestep)

function _concentration_field(pos,model)
    pos_idx = get_spatial_index(pos, model.xmesh, model)
    return model.integrator.u[pos_idx]
end # function
function _concentration_gradient(pos, model)
    pos_idx = get_spatial_index(pos, model.xmesh, model)
    return model.∇u[pos_idx]
end # function
function _concentration_time_derivative(pos, model)
    pos_idx = get_spatial_index(pos, model.xmesh, model)
    return model.∂ₜu[pos_idx]
end # function
                                
model_properties = Dict(
    :xmesh => xs,
    :∇u => ∇u,
    :∂ₜu => ∂ₜu,
    :concentration_field => _concentration_field,
    :concentration_gradient => _concentration_gradient,
    :concentration_time_derivative => _concentration_time_derivative
)

model = initialise_model(;
    microbes, timestep,
    extent, spacing, periodic = false,
    diffeq = true, ode_integrator,
    model_properties
)

@info "initialised model"

adata = [:pos, :state]
u_field(model) = copy(model.integrator.u)
mdata = [u_field]
when = range(0, nsteps; step=round(Int, 5/timestep))
when_model = range(0, nsteps; step=round(Int, 30/timestep))

my_microbe_step!(microbe, model) = microbe_step!(
    microbe, model;
    affect! = brownberg_affect!,
    turnrate = brownberg_turnrate
)

function update_model!(model)
    # update gradient
    finitediff!(model.∇u, model.integrator.u, 1/model.space.spacing)
    # update time derivative 
    get_du!(model.∂ₜu, model.integrator)
end # function

my_model_step!(model) = model_step!(model; update_model!)

adf, mdf = run!(
    model, microbe_step!, my_model_step!, nsteps;
    adata, mdata, when, when_model
)

@info "simulation complete"

linecolors = palette(:plasma, size(mdf,1))
p1 = plot(color_palette = linecolors)
for row in eachrow(mdf)
    plot!(p1, xs[3:end-2], row[:u_field][3:end-2], lw=2, lab=false)
end # for

traj = vectorize_adf_measurement(adf, :pos)
p2 = plot(
    first.(traj)', when.*timestep, lab=false, lw=0.5,
    line_z=when, color=:plasma, colorbar=false
)

plot(p1,p2,xticks=[0,extent/2,extent])