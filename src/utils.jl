export
    rand_vel

"""
    rand_vel(N)
    rand_vel(rng, N)
Generate a random N-tuple of unitary norm.
"""
function rand_vel(D)
    v = rand(Random.default_rng(), D) .* 2 .- 1
    Tuple(v ./ norm(v))
end # function

function rand_vel(rng, D)
    v = rand(rng, D) .* 2 .- 1
    Tuple(v ./ norm(v))
end # function
