export ObstacleSphere, stick!, glide!, bounce!

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
        z = @. x + ε*s
        move_agent!(microbe, z, model)
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
        z = @. x + ε*s
        move_agent!(microbe, z, model)
    end # if
end # function

"""
    bounce!(microbe, sphere::ObstacleSphere, model; ζ=1.0)
`microbe` is reflected off the `sphere` surface, inverting
the direction of its `vel` field.

The parameter `ζ` is the absorption factor of the sphere;
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
        z = @. x + (1+ζ)*ε*s
        move_agent!(microbe, z, model)
        microbe.vel = .-microbe.vel
    end # if
end # function
