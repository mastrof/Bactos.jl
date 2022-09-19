var documenterSearchIndex = {"docs":
[{"location":"#BacteriaBasedModels.jl","page":"Home","title":"BacteriaBasedModels.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Agent-based modelling of bacterial motility and chemotaxis in Julia.","category":"page"},{"location":"#Features","page":"Home","title":"Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Built on Agents.jl\nSimulate bacterial motility in 1, 2 and 3 dimensions with customizable motile patterns.\nSimulate chemotaxis with different tunable models (currently implemented Brown-Berg, Brumley).\nEvaluate quantity of interest (mean-squared displacement, autocorrelation functions...).","category":"page"},{"location":"#Next-steps","page":"Home","title":"Next steps","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Integrate differential equations (concentration fields, flow fields) via DifferentialEquations.jl in parallel with agent stepping.\nImplement more motile/search strategies (e.g. chemokinesis, infotaxis)\nExtend the set of core functionalities (e.g. encounters, interactions)\nImplement complex environments (obstacles, pores)\nInclude models for solid bacteria with arbitrary shapes","category":"page"},{"location":"tutorial/#Tutorial","page":"Tutorial","title":"Tutorial","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Ready-to-run scripts for the functionalities introduced here can be found in the examples directory of the repo.","category":"page"},{"location":"tutorial/#Creating-a-bacterium","page":"Tutorial","title":"Creating a bacterium","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Bacteria are represented by custom types, that must be subtypes of the AbstractAgent type implemented by Agents.jl.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"AbstractMicrobe","category":"page"},{"location":"tutorial/#BacteriaBasedModels.AbstractMicrobe","page":"Tutorial","title":"BacteriaBasedModels.AbstractMicrobe","text":"AbstractMicrobe{D} <: AbstractAgent where {D<:Integer}\nYourMicrobeType{D} <: AbstractMicrobe{D}\n\nAll microbe types in BacteriaBasedModels.jl simulations must be instances of user-defined types that are subtypes of AbstractMicrobe. The parameter D defines the dimensionality of the space in which the microbe type lives (1, 2 and 3 are currently supported).\n\nAll microbe types should have the following fields: - id::Int → id of the microbe - pos::NTuple{D,Float64} → position of the microbe - vel::NTuple{D,Float64} → velocity of the microbe - motility → motile pattern of the microbe - turn_rate::Float64 → average reorientation rate of the microbe - rotational_diffusivity::Float64 → rotational diffusion coefficient\n\n\n\n\n\n","category":"type"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"By default, BacteriaBasedModels provides a basic Microbe type, that is usually sufficient for the simplest types of simulations.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Microbe","category":"page"},{"location":"tutorial/#BacteriaBasedModels.Microbe","page":"Tutorial","title":"BacteriaBasedModels.Microbe","text":"Microbe{D} <: AbstractMicrobe{D}\n\nBasic microbe type for simple simulations.\n\nDefault parameters:\n\nid::Int → identifier used internally by Agents.jl\npos::NTuple{D,Float64} = ntuple(zero,D) → position\nmotility = RunTumble() → motile pattern\nvel::NTuple{D,Float64} = rand_vel(D) .* rand(motility.speed) → velocity vector\nturn_rate::Float64 = 1.0 → frequency of reorientations\nstate::Float64 → generic variable for a scalar internal state\nrotational_diffusivity::Float64 = 0.0 → rotational diffusion coefficient\nradius::Float64 = 0.0 → equivalent spherical radius of the microbe\n\n\n\n\n\n","category":"type"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"In order to create a Microbe living in a one-dimensional space we can just call","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Microbe{1}(id=0)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"It is required to pass a value to the id argument (this behavior might change in the future). All the other parameters will be given default values (as described in the type docstring) if not assigned explicitly.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Similarly, for bacteria living in two or three dimensions we can use","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Microbe{2}(id=0)\nMicrobe{3}(id=0)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Custom parameters can be set via kwargs:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Microbe{3}(\n    id = 0,\n    pos = (300.0, 0.0, 0.0),\n    motility = RunTumble(speed = Normal(40.0, 4.0)),\n    vel = rand_vel(3) .* 40.0,\n    turn_rate = 1.5,\n    state = 0.0,\n    rotational_diffusivity = 0.035,\n    radius = 0.5\n)","category":"page"},{"location":"tutorial/#Creating-a-model","page":"Tutorial","title":"Creating a model","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"BacteriaBasedModel provides a fast way to initialise an AgentBasedModel (from Agents.jl) via the initialise_model function, using a typical procedure. If higher levels of customization are needed, the model will need to be created by hand.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"initialise_model","category":"page"},{"location":"tutorial/#BacteriaBasedModels.initialise_model","page":"Tutorial","title":"BacteriaBasedModels.initialise_model","text":"initialise_model(;\n    microbes,\n    timestep,\n    extent, spacing = extent/20, periodic = true,\n    random_positions = true,\n    model_properties = Dict()\n)\n\nInitialise an AgentBasedModel from population microbes. Requires the integration timestep and the extent of the simulation box.\n\nWhen random_positions = true the positions assigned to microbes are ignored and new ones, extracted randomly in the simulation box, are assigned; if random_positions = false the original positions in microbes are kept.\n\nAny extra property can be assigned to the model via the model_properties dictionary.\n\n\n\n\n\n","category":"function"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"We can now generate a population of microbes and, after choosing an integration timestep and a domain size, we initialise our model, placing the microbes at random locations in the domain.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"microbes = [Microbe{3}(id=i) for i in 1:10]\ntimestep = 0.1\nextent = 100.0\nmodel = initialise_model(;\n    microbes = microbes,\n    timestep = timestep,\n    extent = extent\n)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"AgentBasedModel with 10 agents of type Microbe\n space: periodic continuous space with (100.0, 100.0, 100.0) extent and spacing=5.0\n scheduler: fastest\n properties: timestep","category":"page"},{"location":"tutorial/#Random-walks","page":"Tutorial","title":"Random walks","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Now we can already generate random walks. The setup follows previous sections.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"dt = 0.1\nL = 100.0\nnmicrobes = 8\nmicrobes = [Microbe{1}(id=i, pos=(L/2,)) for i in 1:nmicrobes]\n\nmodel = initialise_model(;\n    microbes = microbes,\n    timestep = dt,\n    extent = L, periodic = false,\n    random_positions = false\n)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Now we need to define the adata variable to choose what observables we want to track, throughout the simulation, for each agent in the system. In our case, only the position field","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"adata = [:pos]","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Now we can run the simulation; the microbe_step! function will take care of the stepping and reorientations:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"nsteps = 1000\nadf, = run!(model, microbe_step!, nsteps; adata)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"x = vectorize_adf_measurement(adf, :pos) .|> first\nplot(\n    (0:nsteps).*dt, x',\n    legend = false,\n    xlab = \"time\",\n    ylab = \"position\"\n)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: One-dimensional random walks of 8 microbes starting from same position)","category":"page"}]
}
