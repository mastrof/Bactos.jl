# BacteriaBasedModels.jl

Agent-based modelling of bacterial motility and chemotaxis in Julia.

## Features
- Built on [Agents.jl](https://juliadynamics.github.io/Agents.jl/stable/)
- Simulate bacterial motility in 1, 2 and 3 dimensions with customizable motile patterns.
- Simulate chemotaxis with different tunable models (currently implemented Brown-Berg, Brumley).
- Evaluate quantity of interest (mean-squared displacement, autocorrelation functions...).

## Next steps
- Integrate differential equations (concentration fields, flow fields) via [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/) in parallel with agent stepping.
- Implement more motile/search strategies (e.g. chemokinesis, infotaxis)
- Extend the set of core functionalities (e.g. encounters, interactions)
- Implement complex environments (obstacles, pores)
- Include models for solid bacteria with arbitrary shapes