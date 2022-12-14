export
    rand_vel, rand_speed, vectorize_adf_measurement

"""
    rand_vel([rng,] N)
Generate a random N-tuple of unitary norm.
"""
function rand_vel(D::Int)
    v = rand(D) .* 2 .- 1
    Tuple(v ./ norm(v))
end # function

function rand_vel(rng, D::Int)
    v = rand(rng, D) .* 2 .- 1
    Tuple(v ./ norm(v))
end # function

"""
    rand_speed(m::AbstractMotilityOneStep)
Extract value from the speed distribution of the motility pattern `m.speed`.
"""
rand_speed(m::AbstractMotilityOneStep) = rand(m.speed)
"""
    rand_speed(m::AbstractMotilityTwoStep)
Extract value from the speed distribution of the motility pattern.
If `motilestate(m) == ForwardState()` extract from `speed_forward`, otherwise
from `speed_backward`.
"""
function rand_speed(m::AbstractMotilityTwoStep)
    s = motilestate(m)
    if s == ForwardState()
        return rand(m.speed_forward)
    else
        return rand(m.speed_backward)
    end
end

"""
    rand_vel([rng,] N::Int, m::AbstractMotility)
Generate a random N-tuple, with norm defined by the speed distribution of `m`.
"""
rand_vel(D::Int, m::AbstractMotility) = rand_vel(D) .* rand_speed(m)
rand_vel(rng, D::Int, m::AbstractMotility) = rand_vel(rng, D) .* rand_speed(m)

"""
    vectorize_adf_measurement(adf, sym)
Collect quantity `sym` from the agent dataframe `adf` and return it in matrix 
form with dimensions (microbes, times).
"""
function vectorize_adf_measurement(adf, sym)
    nmicrobes = unique(adf[!,:id]) |> length
    nsteps = unique(adf[!,:step]) |> length
    datatype = typeof(adf[1,sym])
    s = Matrix{datatype}(undef, nmicrobes, nsteps)
    for t in 1:nsteps
        for i in 1:nmicrobes
            s[i,t] = adf[i + (t-1)*nmicrobes, sym]
        end # for
    end # for
    return s
end # function
