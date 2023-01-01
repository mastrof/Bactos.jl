export
    initialise_model, initialise_ode, initialise_pathfinder

"""
    initialise_model(;
        microbes,
        timestep,
        extent, spacing = extent/20, periodic = true,
        random_positions = true,
        model_properties = Dict(),
        kwargs...
    )
Initialise an `AgentBasedModel` from population `microbes`.
Requires the integration `timestep` and the `extent` of the simulation box.

When `random_positions = true` the positions assigned to `microbes` are
ignored and new ones, extracted randomly in the simulation box, are assigned;
if `random_positions = false` the original positions in `microbes` are kept.

Extra properties can be assigned to the model via the `model_properties` dictionary.

Further kwargs will be sent through to the `AgentBasedModel` constructor from
Agents.jl, so that all its functionalities can be accessed.
"""
function initialise_model(;
    microbes,
    timestep,
    extent, spacing = minimum(extent)/20, periodic = true,
    random_positions = true,
    model_properties = Dict(),
    kwargs...
)
    space_dim = length(microbes[1].pos)
    if typeof(extent) <: Real
        domain = Tuple(fill(extent, space_dim))
    else
        if length(extent) ≠ space_dim
            throw(ArgumentError(
                "Space extent and microbes must have the same dimensionality."
            ))
        end # if
        domain = extent
    end # if

    properties = Dict(
        :t => 0,
        :timestep => timestep,
        :compound_diffusivity => 608.0,
        :concentration_field => (pos,model) -> 0.0,
        :concentration_gradient => (pos,model) -> zero.(pos),
        :concentration_time_derivative => (pos,model) -> 0.0,
        :update! => model_step!,
        model_properties...
    )

    space = ContinuousSpace(domain; spacing, periodic)

    # falls back to eltype(microbes) if there is a single microbe type,
    # builds a Union type if eltype(microbes) is abstract
    MicrobeType = Union{typeof.(microbes)...}

    model = ABM(
        MicrobeType, space;
        properties,
        scheduler = Schedulers.fastest,
        warn = false,
        kwargs...
    )

    for microbe in microbes
        if random_positions
            add_agent!(microbe, model)
        else
            add_agent!(microbe, microbe.pos, model)
        end # if
    end # for

    return model
end # function

"""
    initialise_ode(ode_step!, u₀, p; alg=Tsit5(), kwargs...)
Initialise an OrdinaryDiffEq integrator, using the in-place stepping algorithm
`ode_step!`, initial conditions `u₀` and parameters `p`.
Default integration algorithm is `Tsit5` (others can be accessed by importing
OrdinaryDiffEq).
Any extra parameter can be passed over to the integrator via kwargs.
"""
function initialise_ode(ode_step!, u₀, p; alg=Tsit5(), kwargs...)
    prob = ODEProblem(ode_step!, u₀, (0.0, Inf), p)
    integrator = init(prob, alg; kwargs...)
    return integrator
end # function


function initialise_pathfinder(
    extent::Real, periodic::Bool,
    walkmap::BitArray{D}
) where D
    initialise_pathfinder(ntuple(_->extent,D), periodic, walkmap)
end
function initialise_pathfinder(
    extent::NTuple{D,<:Real}, periodic::Bool,
    walkmap::BitArray{D}
) where D
    space = ContinuousSpace(extent; periodic)
    AStar(space; walkmap)
end