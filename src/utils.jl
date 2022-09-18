export
    rand_vel, vectorize_adf_measurement

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
