export add_encounters!

"""
    add_encounters!(model, bodies;
        key_bodies = :bodies,
        key_encounters = :encounters,
        after_encounter! = reinsert_random!
    )
Add functionalities to determine encounters between the microbes in `model`
and a set of `bodies`.

## Keywords
- `key_bodies = :bodies` assigns a key to `bodies` which will be stored in `model.properties[key_bodies]`
- `key_encounters = :encounters` assigns a key to the number of encounters (`model.properties[key_encounters]`)
- `after_encounter! = reinsert_random!` determines what happens to a microbe after an encounter with a body
"""
function add_encounters!(model::ABM, bodies::AbstractVector;
    key_bodies::Symbol = :bodies,
    key_encounters::Symbol = :encounters,
    after_encounter!::Function = reinsert_random!
)
    if haskey(model.properties, key_bodies)
        @warn "A $key_bodies property already exists; it will be overwritten."
    end
    if haskey(model.properties, key_encounters)
        @warn "A $key_encounters property already exists; it will be overwritten."
    end
    model.properties[key_bodies] = bodies
    model.properties[key_encounters] = 0
    chain!(model,
        (model) -> encounters!(model, key_bodies, key_encounters, after_encounter!)
    )
end

"""
    encounters!(model::ABM, key_bodies, key_encounters, after_encounter!)
Evaluate encounters between all microbes in the `model` and a collection
of bodies stored in `model.properties[key_bodies]`.
Occurrence of an encounter is internally determined by the `is_encounter` function,
with signature `(microbe, body, model)`.
For certain `body` types, `is_encounter` is already defined by Bactos.jl;
if you need extra or different functionality overload `Bactos.is_encounter`.

Whenever an encounter occurs, `model.properties[key_encounters]` is increased by 1.
Then the `after_encounter!` function is applied, which determines what
happens after the encounter (e.g. the microbe can be reinserted or killed).

`encounters!` should be chained to `model.update!`.
"""
function encounters!(model::ABM,
    key_bodies::Symbol, key_encounters::Symbol, after_encounter!::Function
)
    for microbe in allagents(model)
        for body in model.properties[key_bodies]
            if is_encounter(microbe, body, model)
                model.properties[key_encounters] += 1
                after_encounter!(microbe, model, key_bodies)
                break # move to the next microbe
            end
        end
    end
end

"""
    is_encounter(microbe::AbstractMicrobe, sphere::ObstacleSphere, model::ABM)::Bool
Test if an encounter between `microbe` and `sphere` occurred.

## Detailed explanation
More exactly, the function tests whether an encounter will occur at the *next* step.
Since reorientation is the last thing to happen to a microbe at step `t`,
and displacement is the first thing to happen at step `t+1`, at the end of step
`t` we can know exactly whether an encounter will occur during `t+1`.

Instead of only a naive distance check between `microbe` and `sphere` (which would
be unreliable in those cases where `sphere.radius` is comparable or smaller
than the microbial displacement during a single timestep,
`norm(microbe.vel)*model.timestep`), we say that an encounter has occured if
the displacement vector of `microbe` from time `t` to `t+1` intersects
the surface of the `sphere`.

See http://paulbourke.net/geometry/circlesphere/ for geometrical derivation.
"""
function is_encounter(
    microbe::AbstractMicrobe{D}, sphere::ObstacleSphere{D},
    model::ABM
)::Bool where D
    # if microbe is inside sphere, return true 
    is_inside(microbe, sphere) && return true
    # otherwise, compute intersection between displacement and sphere surface
    pos_next = microbe.pos .+ microbe.vel .* model.timestep
    x₁ = microbe.pos .- sphere.pos
    x₂ = pos_next .- sphere.pos
    Δx = x₂ .- x₁
    a = dot(Δx, Δx)
    b = 2.0 * dot(Δx, x₁)
    R = sphere.radius + microbe.radius
    c = dot(x₁, x₁) - R*R
    S = b*b - 4*a*c
    # S<0 → no intersection
    S < 0 && return false
    u₁ = (-b + √S) / (2a)
    u₂ = (-b - √S) / (2a)
    # the segment intersects the sphere if at least one u ∈ [0,1]
    return (0 ≤ u₁ ≤ 1) || (0 ≤ u₂ ≤ 1)
end


"""
    reinsert_random!(microbe::AbstractMicrobe, model::ABM, bodies::Symbol; avoid::Bool = true)
Reinsert `microbe` in a random position.
If `avoid` is set to `true`, only positions which don't trigger
a new encounter with any body in `model.properties[bodies]` are accepted.
"""
function reinsert_random!(
    microbe::AbstractMicrobe, model::ABM, key_bodies::Symbol;
    avoid::Bool = true
)
    while true
        any_overlap = false
        # move to a new position at random
        move_agent!(microbe, model)
        # don't perform any check if avoid==false
        (~avoid) && break
        # if avoid==true, extract new position until it triggers no encounter
        for body in model.properties[key_bodies]
            if is_encounter(microbe, body, model)
                # break if position is not valid
                any_overlap = true
                break
            end
            any_overlap && break # extract a new position if there is an overlap
        end
        # break out of while loop if position is valid
        (~any_overlap) && break
    end
end
