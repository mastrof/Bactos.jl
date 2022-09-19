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
