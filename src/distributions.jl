export Degenerate, Arccos

struct Degenerate{T<:Real} <: ContinuousUnivariateDistribution
    x::T
end # struct
Base.rand(d::Degenerate) = d.x

struct Arccos{T<:Real} <: ContinuousUnivariateDistribution
    a::T
    b::T
    Arccos{T}(a::T, b::T) where {T} = new{T}(a::T, b::T)
end
function Arccos(a::Real, b::Real; check_args::Bool=true)
    Distributions.@check_args Arccos a<b -1≤a≤1 -1≤b≤1
    return Arccos{Float64}(Float64(a), Float64(b))
end # function
Arccos() = Arccos(-1,1)
Base.rand(rng::AbstractRNG, d::Arccos) = acos(rand(rng, Uniform(d.a, d.b)))