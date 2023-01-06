export
    add_neighborlist!, init_neighborlist, update_neighborlist!,
    surface_interaction!


function add_neighborlist!(model::ABM,
    x::AbstractVector, cutoff::Real;
    listname::Symbol = :neighborlist, x_or_y::Char = 'x'
)
    neighborlist = init_neighborlist(x, model.space.extent, cutoff, model.space.periodic)
    model.properties[listname] = neighborlist
    chain!(model, m::ABM -> update_neighborlist!(m; listname, x_or_y))
end
function add_neighborlist!(model::ABM,
    x::AbstractVector, y::AbstractVector, cutoff::Real;
    listname::Symbol = :neighborlist, x_or_y::Char = 'x'
)
    neighborlist = init_neighborlist(x, y, model.space.extent, cutoff, model.space.periodic)
    model.properties[listname] = neighborlist
    chain!(model, m::ABM -> update_neighborlist!(m; listname, x_or_y))
end

"""
init_neighborlist(
    x::AbstractVector,
    extent::Union{<:Real,NTuple{D,<:Real}},
    cutoff::Real, periodic::Bool
) where D
Initialise a neighbor list between objects in `x`
in a domain of size `extent` with a neighbor cutoff radius given by `cutoff`.
`periodic` determines whether the system has periodic boundary conditions.
Uses the `PeriodicSystem` interface of `CellListMap`.
"""
function init_neighborlist(
    x::AbstractVector,
    extent::Real, cutoff::Real, periodic::Bool
)
    D = length(getpos(first(x)))
    init_neighborlist(x, ntuple(_->extent,D), cutoff, periodic)
end
function init_neighborlist(
    x::AbstractVector,
    extent::NTuple{D,<:Real}, cutoff::Real, periodic::Bool
) where D
    xpos = getpos.(x)
    return PeriodicSystem(
        xpositions = xpos,
        unitcell = SVector(extent) .+ (periodic ? 0.0 : cutoff),
        cutoff = cutoff,
        output = 0.0
    )
end # function

"""
    init_neighborlist(
        x::AbstractVector, y::AbstractVector,
        extent::Union{<:Real,NTuple{D,<:Real}},
        cutoff::Real, periodic::Bool
    ) where D
Initialise a neighbor list between objects in `x` and `y`
in a domain of size `extent` with a neighbor cutoff radius given by `cutoff`.
`periodic` determines whether the system has periodic boundary conditions.
Uses the `PeriodicSystem` interface of `CellListMap`.
"""
function init_neighborlist(
    x::AbstractVector, y::AbstractVector,
    extent::Real, cutoff::Real, periodic::Bool
)
    D = length(getpos(first(x)))
    init_neighborlist(x, y, ntuple(_->extent,D), cutoff, periodic)
end
function init_neighborlist(
    x::AbstractVector, y::AbstractVector,
    extent::NTuple{D,<:Real}, cutoff::Real, periodic::Bool
) where D
    xpos = getpos.(x)
    ypos = getpos.(y)
    return PeriodicSystem(
        xpositions = xpos,
        ypositions = ypos,
        unitcell = SVector(extent) .+ (periodic ? 0.0 : cutoff),
        cutoff = cutoff,
        output = 0.0
    )
end # function

getpos(x::AbstractVector) = SVector{length(x)}(x)
getpos(x::NTuple{D,<:Real}) where D = SVector{D}(x)
getpos(x::AbstractMicrobe{D}) where D = SVector{D}(x.pos)
getpos(x::ObstacleSphere{D}) where D = SVector{D}(x.pos)


"""
    update_neighborlist!(model; listname=:neighborlist, x_or_y='x')
Update the position of all microbes in a neighbor list in `model`.
Must be passed in the `update_model!` function to `model_step!`.

## Keywords
* `listname::Symbol = :neighborlist`: name of the neighbor list to update
    (`model.properties[listname]`); this is useful if different neighbor lists are used
    (e.g. to evaluate different types of interactions).
* `x_or_y::Char = 'x'`: update the `xpositions` or the `ypositions` field in
    the neighbor list.
"""
function update_neighborlist!(model::ABM;
    listname::Symbol = :neighborlist,
    x_or_y::Char = 'x'
)
    for microbe in allagents(model)
        update_neighborlist!(microbe, model; listname, x_or_y)
    end
end
"""
    update_neighborlist!(microbe, model; listname=:neighborlist, x_or_y='x')
Update the position of `microbe` in a neighbor list.
Must be passed as the `affect!` function to `microbe_step!`.

## Keywords
* `listname::Symbol = :neighborlist`: name of the neighbor list to update
    (`model.properties[listname]`); this is useful if different neighbor lists are used
    (e.g. to evaluate different types of interactions).
* `x_or_y::Char = 'x'`: update the `xpositions` or the `ypositions` field in
    the neighbor list.
"""
function update_neighborlist!(microbe::AbstractMicrobe, model::ABM;
    listname::Symbol = :neighborlist,
    x_or_y::Char = 'x'
)
    neighborlist = model.properties[listname]
    if lowercase(x_or_y) == 'x'
        neighborlist.xpositions[microbe.id] = SVector(microbe.pos)
    elseif lowercase(x_or_y) == 'y'
        neighborlist.ypositions[microbe.id] = SVector(microbe.pos)
    else
        throw(ArgumentError(
            "Value $(x_or_y) not valid for `x_or_y` keyword."
        ))
    end
    return nothing
end # function


function surface_interaction!(x,y,i,j,d²,f,model;bodies=:bodies)
    body = model.property[bodies][j]
    body.affect!(model[i], body, model)
    return f
end # function

"""
    surface_interaction!(model; listname=:neighborlist, bodies=:bodies)
Evaluate the effect of surface interactions between microbes and bodies
using the neighborlist for efficient computation.
Requires the neighborlist (initialised via `init_neighborlist`) to be set
as a model property `model.properties[listname]`,
and a list of bodies with which interactions occur (`model.properties[bodies]`).
"""
function surface_interaction!(model::ABM;
    listname::Symbol=:neighborlist, bodies::Symbol=:bodies
)
    map_pairwise!(
        (x,y,i,j,d²,f) -> surface_interaction!(x,y,i,j,d²,f,model;bodies),
        model.properties[listname]
    )
    return nothing
end # function