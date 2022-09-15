export
    rand_vel

function rand_vel(D)
    v = rand(Random.default_rng(), D) .* 2 .- 1
    Tuple(v ./ norm(v))
end # function

function rand_vel(rng, D)
    v = rand(rng, D) .* 2 .- 1
    Tuple(v ./ norm(v))
end # function