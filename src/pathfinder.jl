export add_pathfinder!, initialise_pathfinder, pathfinder_step!

"""
    add_pathfinder!(model::ABM, walkmap::BitArray)
Add pathfinding to `model` in the domain described by `walkmap`.
`walkmap` must have the same dimensions as `model.space`
(but its resolution is unrelated to `model.space.spacing`)
"""
function add_pathfinder!(model::ABM, walkmap::BitArray)
    pathfinder = initialise_pathfinder(model.space, walkmap)
    model.properties[:pathfinder] = pathfinder
end

function initialise_pathfinder(
    extent::Real, periodic::Bool,
    walkmap::BitArray{D}
) where D
    initialise_pathfinder(ntuple(_->extent,D), periodic, walkmap)
end
function initialise_pathfinder(
    extent::NTuple{D,<:Real}, periodic::Bool,
    walkmap::BitArray{D}
) where D
    space = ContinuousSpace(extent; periodic)
    AStar(space; walkmap)
end
function initialise_pathfinder(space::ContinuousSpace{D}, walkmap::BitArray{D}) where D
    AStar(space; walkmap)
end

"""
    pathfinder_step!(microbe::AbstractMicrobe, model::ABM, dt::Real)
Perform an integration step for `microbe` motion with pathfinding
(`model.pathfinder`).
"""
function pathfinder_step!(microbe::AbstractMicrobe, model::ABM, dt::Real)
    target_position = microbe.pos .+ microbe.vel .* dt
    U = norm(microbe.vel)
    plan_route!(microbe, target_position, model.pathfinder)
    move_along_route!(microbe, model, model.pathfinder, U, dt)
end
