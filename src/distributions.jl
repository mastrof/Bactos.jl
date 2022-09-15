export Degenerate

struct Degenerate{T<:Real} <: ContinuousUnivariateDistribution
    x::T
end # struct

Base.rand(d::Degenerate) = d.x