export model_step!

"""
    model_step!(model; update_model! = (model) -> nothing)
Update model properties through the `update_model!` function
(defaults to nothing).
If `model` contains an OrdinaryDiffEq integrator among its
properties (`model.integrator`), also perform an integration step.
"""
function model_step!(model;
    update_model!::Function = (model) -> nothing
)
    # if a diffeq integrator is provided, integrate over a timestep
    if haskey(model.properties, :integrator)
        step!(model.integrator, model.timestep, true)
    end # if
    # increase step count
    model.t += 1
    # update model properties
    update_model!(model)
    return nothing
end # function