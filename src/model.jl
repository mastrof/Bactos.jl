export initialise_model, model_step!

# convenience construct to dispatch on ABM with microbe agents
BBM = ABM{
    <:ContinuousSpace,
    <:AbstractMicrobe,
    <:Function,
    <:AbstractDict,
    <:AbstractRNG
}

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
    model_step!(model::ABM)
Model stepping function, stored in `model.update!`.
By default, it only increases the time count of `model` (`model.t += 1`).

Extra functionalities can be added to the basic `model_step!` through
function chaining e.g. `chain!(model, (f₁!, f₂!, f₃!))`.

Some add-ons of general utility are already provided by Bactos.jl
(see respective documentations for more details):
* `diff_step!`: used to integrate differential equations if `model` has an
    `integrator` property (`model.integrator`)
* `update_neighborlist!`: update the microbe positions in a neighborlist
* `surface_interaction!`: used to evaluate the effect of interactions between
    microbes and other surfaces
"""
function model_step!(model::ABM)
    model.t += 1
    return nothing
end # function
