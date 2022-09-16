export
    microbe_step!

function microbe_step!(
    microbe, model;
    affect! = (microbe, model) -> nothing,
    turnrate = (microbe, model) -> microbe.turn_rate
)
    dt = model.timestep # integration timestep
    # update microbe position
    move_agent!(microbe, model, dt)
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

