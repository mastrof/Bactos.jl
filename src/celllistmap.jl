export
    init_neighborlist, update_microbepos_neighborlist!,
    surface_interaction!

"""
    init_neighborlist(microbes::Vector, bodies::Vector, extent, cutoff::Float64, periodic::Bool)
Initialise a neighbor list using the `PeriodicSystem` interface of
`CellListMap`, between `x` and `y` in a domain of size `extent` with
a neighbor cutoff radius given by `cutoff`.
"""
function init_neighborlist(microbes, bodies, extent, periodic, cutoff)
    microbepos = [SVector(m.pos) for m in microbes]
    bodypos = [SVector(b.pos) for b in bodies]
    return PeriodicSystem(
        xpositions = microbepos,
        ypositions = bodypos,
        unitcell = SVector(extent) .+ (periodic ? 0.0 : cutoff),
        cutoff = cutoff,
        output = 0.0
    )
end # function

"""
    update_microbepos_neighborlist!(microbe, model)
Update the position of `microbe` in `model.neighborlist`;
required for correct evaluation of interactions.
This function must be passed as the `affect!` function to
`microbe_step!`.
"""
function update_microbepos_neighborlist!(microbe, model)
    model.neighborlist.xpositions[microbe.id] = SVector(microbe.pos)
end # function

function surface_interaction!(x,y,i,j,d²,f,model)
    body = model.bodies[j]
    body.affect!(model[i], body, model)
    return f
end # function

"""
    surface_interaction!(model)
Evaluate the effect of surface interactions between microbes and bodies
using the neighborlist for efficient computation.
Requires the neighborlist (initialised via `init_neighborlist`) to be set
as a model property `model.neighborlist`.
"""
function surface_interaction!(model)
    map_pairwise!(
        (x,y,i,j,d²,f) -> surface_interaction!(x,y,i,j,d²,f,model),
        model.neighborlist
    )
end # function