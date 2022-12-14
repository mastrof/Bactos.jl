export 
    rotate, turn!, rotational_diffusion!

rotate(w::SVector{1}) = -w
rotate(w::SVector{1}, θ) = rotate(w)
rotate(w::SVector{1}, θ, ϕ) = rotate(w)

rotate(w::SVector{2}, θ) = Angle2d(θ) * w
rotate(w::SVector{2}, θ, ϕ) = rotate(w, θ)

function rotate(w::SVector{3}, θ, ϕ)
    m = findfirst(w .≠ 0)
    n = m%3 + 1
    u = SVector(0., 0., 0.)
    u = setindex(u, w[m], n)
    u = setindex(u, -w[n], m)
    # rotate w around its normal u 
    a = AngleAxis(θ, u...) * w
    # rotate a around the original w direction
    return AngleAxis(ϕ, w...) * a
end # function 

rotate(w::Tuple, θ, ϕ) = rotate(SVector(w), θ, ϕ)
rotate(w::Tuple, θ) = rotate(SVector(w), θ)


function turn!(microbe::AbstractMicrobe, motility::AbstractMotilityOneStep)
    # store actual speed
    U₀ = norm(microbe.vel)
    # perform reorientation
    θ = rand(motility.polar)
    ϕ = rand(motility.azimuthal)
    microbe.vel = rotate(microbe.vel, θ, ϕ) |> Tuple
    # extract new speed from distribution
    U₁ = rand(motility.speed)
    # update speed
    microbe.vel = microbe.vel .* (U₁ / U₀)
    return nothing
end # function

function turn!(microbe::AbstractMicrobe, motility::AbstractMotilityTwoStep)
    # store current speed
    U₀ = norm(microbe.vel)
    # perform reorientation depending on current motile state
    if motilestate(motility) == ForwardState()
        # reorient according to forward-state angles
        θ = rand(motility.polar_forward)
        ϕ = rand(motility.azimuthal_forward)
        # extract new speed from backward-state distribution
        U₁ = rand(motility.speed_backward)
    elseif motilestate(motility) == BackwardState()
        # reorient according to backward-state angles
        θ = rand(motility.polar_backward)
        ϕ = rand(motility.azimuthal_backward)
        # extract new speed from forward-state distribution
        U₁ = rand(motility.speed_forward)
    end # if
    # reorient
    microbe.vel = rotate(microbe.vel, θ, ϕ) |> Tuple
    # update motile state
    switch!(motility.motile_state)
    # update speed
    microbe.vel = microbe.vel .* (U₁ / U₀)
    return nothing
end # function


rotational_diffusion!(microbe::AbstractMicrobe{1}, dt) = nothing
function rotational_diffusion!(microbe::AbstractMicrobe, dt)
    D_rot = microbe.rotational_diffusivity
    σ = sqrt(2*D_rot*dt)
    polar = rand(Normal(0, σ))
    azimuthal = rand(Arccos())
    microbe.vel = Tuple(rotate(microbe.vel, polar, azimuthal))
    return nothing
end # function 