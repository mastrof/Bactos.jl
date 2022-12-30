export ObstacleSphere, get_walkmap, stick!, glide!, bounce!, is_encounter

"""
    struct ObstacleSphere{D}
        pos::NTuple{D,Float64}
        radius::Float64
        affect!::Function
    end

Spherical obstacle in a D-dimensional space.

- `pos`: center position of the sphere
- `radius`: radius of the sphere
- `affect!`: function with signature `(microbe,sphere,model)`, defines microbe-sphere interaction.

Three `affect!` functions are provided out-of-the box:
`stick!`, `glide!` and `bounce!`.
"""
struct ObstacleSphere{D}
    pos::NTuple{D,Float64}
    radius::Float64
    affect!::Function
end # struct
ObstacleSphere(pos::NTuple{D,<:Real}, radius::Real,
    affect!::Function = (_,_,_) -> nothing
) where D = ObstacleSphere{D}(Float64.(pos), Float64(radius), affect!)


function initialise_pathfinder(
    extent, periodic::Bool,
    r::Real, spheres::AbstractVector{ObstacleSphere{D}};
    Δ::Real=r/2
) where D
    walkmap = get_walkmap(extent, r, spheres; Δ)
    initialise_pathfinder(extent, periodic, walkmap)
end

"""
    is_inside(m::AbstractMicrobe, s::ObstacleSphere)::Bool
Check if microbe `m` is inside the volume of sphere `s`.
"""
is_inside(m::AbstractMicrobe, s::ObstacleSphere)::Bool = is_inside(m.pos, m.radius, s)
"""
    is_inside(pos::NTuple{D,Float64}, r::Real, s::ObstacleSphere{D})::Bool where D
Check if a sphere of radius `r` centered at `pos` is inside the volume of sphere `s`.
"""
function is_inside(pos::NTuple{D,Float64}, r::Real, s::ObstacleSphere{D})::Bool where D
    norm(pos .- s.pos) < r + s.radius
end


function get_walkmap(
    extent::NTuple{D,<:Real}, r::Real,
    spheres::AbstractVector{ObstacleSphere{D}};
    Δ::Real=r/2
)::BitArray{D} where D
    mesh = ntuple(i -> 0:Δ:extent[i], D)
    itr = Iterators.product(mesh...)
    return BitArray([is_walkable(pos, r, spheres) for pos in itr])
end
function is_walkable(
    pos::NTuple{D,<:Real}, r::Real,
    spheres::AbstractVector{ObstacleSphere{D}}
)::Bool where D
    for sphere in spheres
        if is_inside(pos, r, sphere)
            return false
        end
    end
    return true
end


"""
    stick!(microbe, sphere::ObstacleSphere, model)
`microbe` sticks to the `sphere` surface at the point of contact.
The orientation of the microbe is unchanged.
"""
function stick!(microbe, sphere::ObstacleSphere, model)
    x = microbe.pos
    y = sphere.pos
    R = microbe.radius + sphere.radius
    d = norm(x .- y)
    if d < R
        s = @. -microbe.vel * model.timestep
        a = sum(abs2.(s))
        b = 2.0 * dot(x.-y, s)
        c = d*d - R*R
        ε = -b/2a * (1 - sqrt(1 - 4*a*c/(b*b)))
        z = @. ε*s
        walk!(microbe, z, model)
    end # if
end # function

"""
    glide!(microbe, sphere::ObstacleSphere, model)
`microbe` sticks to the `sphere` surface while gliding along.
"""
function glide!(microbe, sphere::ObstacleSphere, model)
    x = microbe.pos
    y = sphere.pos
    R = microbe.radius + sphere.radius
    d = norm(x .- y)
    if d < R
        s = y .- x
        a = sum(abs2.(s))
        c = d*d - R*R
        ε = 1 - sqrt(1 - c/a)
        z = @. ε*s
        walk!(microbe, z, model)
    end # if
end # function

"""
    bounce!(microbe, sphere::ObstacleSphere, model; ζ=1.0)
`microbe` is reflected off the `sphere` surface, inverting
the direction of its `vel` field.
The parameter `ζ` is the elastic coefficient of the collision;
for `ζ=1` the collision is perfectly elastic (microbe run length is conserved);
for `ζ=0` the microbe sticks to the surface (but vel is inverted).
"""
function bounce!(microbe, sphere::ObstacleSphere, model; ζ=1.0)
    if !(0 ≤ ζ ≤ 1) 
        throw(DomainError(ζ, "ζ must have value in the range [0,1]."))
    end # if
    x = microbe.pos
    y = sphere.pos
    R = microbe.radius + sphere.radius
    d = norm(x .- y)
    if d < R
        s = @. -microbe.vel * model.timestep
        a = sum(abs2.(s))
        b = 2.0 * dot(x.-y, s)
        c = d*d - R*R
        ε = -b/2a * (1 - sqrt(1 - 4*a*c/(b*b)))
        z = @. (1+ζ)*ε*s
        # hit surface
        z₁ = @. ε*s
        walk!(microbe, z₁, model)
        # reorient
        r_hat = (y.-x)./d
        deflection = -2 .* dot(microbe.vel, r_hat) .* r_hat
        microbe.vel = @. ζ * (microbe.vel + deflection)
        # bounce
        z₂ = @. ε * (microbe.vel * model.timestep)
        walk!(microbe, z₂, model)
    end # if
end # function


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
the surface of the `sphere`).
See http://paulbourke.net/geometry/circlesphere/ for geometrical derivation.
"""
function is_encounter(
    microbe::AbstractMicrobe, sphere::ObstacleSphere,
    model::ABM
)::Bool
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