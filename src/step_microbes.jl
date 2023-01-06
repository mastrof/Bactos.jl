export simulate!, microbe_step!, turnrate, affect!

# does not require microbe_step! and model.update! to be specified
function Agents.run!(model::BBM,
    n = 1;
    when = true,
    when_model = when,
    adata = nothing,
    mdata = nothing,
    obtainer = identity,
    agents_first = true,
    showprogress = false,
)
    run!(model, microbe_step!, model.update!, n;
        when, when_model,
        adata, mdata,
        obtainer, agents_first, showprogress
    )
end

function microbe_step!(microbe::AbstractMicrobe{D}, model::ABM) where D
    dt = model.timestep # integration timestep
    # update microbe position
    if haskey(model.properties, :pathfinder)
        pathfinder_step!(microbe, model, dt)
    else
        move_agent!(microbe, model, dt)
    end
    # reorient microbe due to rotational diffusion
    rotational_diffusion!(microbe, dt)
    # update microbe state
    affect!(microbe, model)
    # update reorientation rate
    ω = turnrate(microbe, model)
    if rand() < ω*dt # if true reorient microbe
        turn!(microbe, microbe.motility)
    end # if
    return nothing
end # function

turnrate(microbe::AbstractMicrobe, model) = microbe.turn_rate
affect!(microbe::AbstractMicrobe, model) = nothing