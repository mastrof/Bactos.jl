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
    u = zeros(3)
    u[n] = w[m]
    u[m] = -w[n]
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
    # store actual speed
    U₀ = norm(microbe.vel)
    # perform reorientation depending on current motile state
    if motility.motile_state[1] == 0
        θ = rand(motility.polar_0)
        ϕ = rand(motility.azimuthal_0)
    else
        θ = rand(motility.polar_1)
        ϕ = rand(motility.azimuthal_1)
    end # if
    microbe.vel = rotate(microbe.vel, θ, ϕ) |> Tuple
    # update motile state
    motility.motile_state[1] = 1 - motility.motile_state[1]
    # extract new speed from distribution
    U₁ = rand(motility.speed)
    # update speed
    microbe.vel = microbe.vel .* (U₁ / U₀)
    return nothing
end # function


rotational_diffusion!(microbe::AbstractMicrobe{1}, dt) = nothing
function rotational_diffusion!(microbe::AbstractMicrobe, dt)
    D_rot = microbe.rotational_diffusivity
    σ = sqrt(2*D_rot*dt)
    #== this might be wrong ==#
    polar = rand(Normal(0, σ))
    azimuthal = rand(Uniform(0, 2π))
    #==#
    microbe.vel = Tuple(rotate(microbe.vel, polar, azimuthal))
    return nothing
end # function 