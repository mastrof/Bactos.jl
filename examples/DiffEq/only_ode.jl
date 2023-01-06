using Bactos
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

concentration_field(x,C,σ,x₀) = C * exp(-(x-x₀)^2/(2*σ^2))
concentration_gradient(x,C,σ,x₀) = -(x-x₀)/σ^2 * concentration_field(x,C,σ,x₀)

nmicrobes = 10
microbes = [Microbe{1}(id=i) for i in 1:nmicrobes]

extent = 100.0
spacing = 0.2
xs = range(-2*spacing, extent+2*spacing; step=spacing)
timestep = 0.1

x₀ = extent/2
C = 5.0
σ = 3.0
concentration_field(x) = concentration_field(x,C,σ,x₀)
concentration_gradient(x) = concentration_gradient(x,C,σ,x₀)
u₀ = concentration_field.(xs)

function odestep!(du, u, p, t)
    D, _dx, = p
    a = D * _dx * _dx
    laplacian!(du, u, a)
    # absorbing walls
    du[1] = du[2] = du[end] = du[end-1] = 0.0
end # function

model = initialise_model(;
    microbes, timestep,
    extent, spacing,
)
add_diffeq!(model, odestep!, u₀, (1.0, 1/spacing); dtmax=spacing^2/2)

u_field(model) = copy(model.integrator.u)
mdata = [u_field]
nsteps = round(Int, 600 / timestep)
when_model = range(0, nsteps; step=round(Int, 30/timestep))
_, mdf = run!(model, nsteps; mdata, when_model)

linecolors = palette(:plasma, size(mdf,1))
plot(color_palette = linecolors)
for row in eachrow(mdf)
    plot!(xs[3:end-2], row[:u_field][3:end-2], lw=2, lab=false)
end # for
plot!()