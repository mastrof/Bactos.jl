export
    rand_vel, rand_speed,
    chain!, chain,
    vectorize_adf_measurement

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
    if m.state == Forward
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
    chain!(model::ABM, gs)
Chain a sequence of functions `gs` to the `model.update!` function.
Each function in `gs` should only take `model` as argument, and 
preferably return `nothing`.

These functions will be applied in the specified order (left-to-right),
*after* the `model.t` field has been increased.
"""
chain!(model::ABM, gs) = model.update! = chain(model.update!, gs)
"""
    chain(f::Function, gs)
Chain function `f` to a collection of functions `gs`
`f` and all the functions in `gs` must accept the same set of arguments.
"""
function chain(f::Function, gs)
    js = eachindex(gs)
    if length(js) > 1
        return chain(f, chain(first(gs), gs[js[2:end]]))
    else
        return chain(f, first(gs))
    end
end
"""
    chain(f::Function, g::Function)
Return a function that successively applies functions `f` and `g`.
`f` and `g` must accept the same set of arguments.
"""
chain(f::Function, g::Function) = (a...;kw...) -> (f(a...;kw...); g(a...;kw...))

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
