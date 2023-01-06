export add_diffeq!, initialise_ode, diffeq_step!

"""
    add_diffeq!(model::ABM, f, u₀, p; alg=Tsit5(), kwargs...)
Add a differential equation integrator to `model`.

## Arguments
* `f`: function that defines the differential equation `du/dt = f`
* `u₀`: initial conditions
* `p`: extra parameters accessible by `f`

## Keywords
* `alg = Tsit5()`: solver algorithm for the differential equation;
    `Tsit5` is used as a default choice, others can be accessed by
    importing OrdinaryDiffEq.jl
* Other keyword arguments are directly passed to the integrator.
"""
function add_diffeq!(model::ABM, f, u₀, p; alg=Tsit5(), kwargs...)
    integrator = initialise_ode(f, u₀, p; alg, kwargs...)
    model.properties[:integrator] = integrator
    chain!(model, diffeq_step!)
    return nothing
end

"""
    initialise_ode(f, u₀, p; alg=Tsit5(), kwargs...)
Initialise an OrdinaryDiffEq integrator, using the function
`f`, initial conditions `u₀` and parameters `p`.
Default integration algorithm is `Tsit5` (others can be accessed by importing
OrdinaryDiffEq).
Any extra parameter can be passed over to the integrator via kwargs.
"""
function initialise_ode(f, u₀, p; alg=Tsit5(), kwargs...)
    prob = ODEProblem(f, u₀, (0.0, Inf), p)
    integrator = init(prob, alg; kwargs...)
    return integrator
end # function

"""
    diffeq_step!(model::ABM)
Integrate a differential equation over a timestep (`model.timestep`).
The differential equation has to be defined through the integrator interface
of DifferentialEquations.jl and stored in `model.integrator`.
"""
diffeq_step!(model::ABM) = step!(model.integrator, model.timestep, true)
