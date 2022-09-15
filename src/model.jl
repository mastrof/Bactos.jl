export
    initialise_model

function initialise_model(;
    microbes,
    timestep,
    extent, spacing = extent/20, periodic = true,
    random_positions = true,
    model_properties = Dict(),
)
    properties = Dict(
        :timestep => timestep,
        model_properties...
    )

    space_dim = length(microbes[1].pos)
    space = ContinuousSpace(
        Tuple(fill(extent, space_dim)),
        spacing = spacing,
        periodic = periodic
    )

    MicrobeType = eltype(microbes)

    model = ABM(
        MicrobeType, space;
        properties,
        scheduler = Schedulers.fastest,
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