export ObstacleSphere, get_walkmap

"""
    struct ObstacleSphere{D}
        pos::NTuple{D,Float64}
        radius::Float64
    end

Spherical obstacle in a D-dimensional space.

- `pos`: center position of the sphere
- `radius`: radius of the sphere
"""
struct ObstacleSphere{D}
    pos::NTuple{D,Float64}
    radius::Float64
end # struct
ObstacleSphere(pos::NTuple{D,<:Real}, radius::Real) where D =
    ObstacleSphere{D}(Float64.(pos), Float64(radius))

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
        if !is_walkable(pos, r, sphere)
            return false
        end
    end
    return true
end
function is_walkable(
    pos::NTuple{D,<:Real}, r::Real,
    sphere::ObstacleSphere{D}
)::Bool where D
    norm(pos .- sphere.pos) ≥ r + sphere.radius
end
