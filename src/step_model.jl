export model_step!, diffeq_step!

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

"""
    diffeq_step!(model::ABM)
Integrate a differential equation over a timestep (`model.timestep`).
The differential equation has to be defined through the integrator interface
of DifferentialEquations.jl and stored in `model.integrator`.
"""
diffeq_step!(model::ABM) = step!(model.integrator, model.timestep, true)