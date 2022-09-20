# Tutorial
Ready-to-run scripts for the functionalities introduced here can be found in the `examples` directory of the [repo](https://github.com/mastrof/BacteriaBasedModels).

## Creating a bacterium
Bacteria are represented by custom types, that must be subtypes of the `AbstractAgent` type implemented by Agents.jl.
```@docs
AbstractMicrobe
```
By default, BacteriaBasedModels provides a basic `Microbe` type, that is usually sufficient for the simplest types of simulations.
```@docs
Microbe
```

In order to create a `Microbe` living in a one-dimensional space we can just call
```julia
Microbe{1}(id=0)
```
It is *required* to pass a value to the `id` argument (this behavior might change in the future).
All the other parameters will be given default values (as described in the type docstring) if not assigned explicitly.

Similarly, for bacteria living in two or three dimensions we can use
```julia
Microbe{2}(id=0)
Microbe{3}(id=0)
```

Custom parameters can be set via kwargs:
```julia
Microbe{3}(
    id = 0,
    pos = (300.0, 0.0, 0.0),
    motility = RunTumble(speed = Normal(40.0, 4.0)),
    vel = rand_vel(3) .* 40.0,
    turn_rate = 1.5,
    state = 0.0,
    rotational_diffusivity = 0.035,
    radius = 0.5
)
```

## Creating a model
BacteriaBasedModel provides a fast way to initialise an `AgentBasedModel` (from Agents.jl) via the
`initialise_model` function, using a typical procedure.
If higher levels of customization are needed, the model will need to be created by hand.
```@docs
initialise_model
```

We can now generate a population of microbes and, after choosing an integration timestep and a domain size, we initialise our model, placing the microbes at random locations in the domain.
```julia
microbes = [Microbe{3}(id=i) for i in 1:10]
timestep = 0.1
extent = 100.0
model = initialise_model(;
    microbes = microbes,
    timestep = timestep,
    extent = extent
)
```

```
AgentBasedModel with 10 agents of type Microbe
 space: periodic continuous space with (100.0, 100.0, 100.0) extent and spacing=5.0
 scheduler: fastest
 properties: timestep
```


## Random walks
Now we can already generate random walks.
The setup follows previous sections.
```julia
dt = 0.1
L = 100.0
nmicrobes = 8
microbes = [Microbe{1}(id=i, pos=(L/2,)) for i in 1:nmicrobes]

model = initialise_model(;
    microbes = microbes,
    timestep = dt,
    extent = L, periodic = false,
    random_positions = false
)
```

Now we need to define the `adata` variable to choose what observables we want to track, throughout the simulation, for each agent in the system. In our case, only the position field
```julia
adata = [:pos]
```
Now we can run the simulation; the `microbe_step!` function will take care of the stepping and reorientations:
```julia
nsteps = 1000
adf, = run!(model, microbe_step!, nsteps; adata)
```

```julia
x = vectorize_adf_measurement(adf, :pos) .|> first
plot(
    (0:nsteps).*dt, x',
    legend = false,
    xlab = "time",
    ylab = "position"
)
```
![One-dimensional random walks of 8 microbes starting from same position](random_walk_1d.png)

Similarly for a two-dimensional random walk, using run-reverse-flick motility and non-zero rotational diffusion:
```julia
dt = 0.1
L = 100.0
nmicrobes = 1
microbes = [
    Microbe{2}(
        id=i, pos=(L/2,L/2),
        motility=RunReverseFlick(),
        rotational_diffusivity = 0.02,
        ) for i in 1:nmicrobes
]

model = initialise_model(;
    microbes = microbes,
    timestep = dt,
    extent = extent, periodic = false,
    random_positions = false,
)

nsteps = 500
adata = [:pos]
adf, = run!(model, microbe_step!, nsteps; adata)

traj = vectorize_adf_measurement(adf, :pos)
x = first.(traj)
y = last.(traj)
plot(
    x', y', line_z = (0:nsteps).*dt,
    legend=false,
    xlab = "x", ylab = "y",
    colorbar = true, colorbar_title = "time"
)
```
![Two-dimensional random walk using run-reverse-flick motility](random_walk_2d.png)


Microbes with different motile patterns can also be combined in the same simulation, without extra complications or computational costs:
```julia
n = 3
microbes_runtumble = [Microbe{2}(id=i, motility=RunTumble()) for i in 1:n]
microbes_runrev = [Microbe{2}(id=n+i, motility=RunReverse()) for i in 1:n]
microbes_runrevflick = [Microbe{2}(id=2n+1, motility=RunReverseFlick()) for i in 1:n]
microbes = vcat(
    microbes_runtumble, microbes_runrev, microbes_runrevflick
)
```