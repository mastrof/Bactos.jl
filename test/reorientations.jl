using Test, Bactos, Random
using LinearAlgebra: norm, dot
using StaticArrays
using Distributions: Uniform, Normal

@testset "Microbe reorientations" begin
    ≊(x,y) = isapprox(x, y, atol=1e-12) # \approxeq

    @testset "One-dimensional reorientations" begin
        # all reorientations in 1D should be just reversals
        motility = RunTumble(speed=[1.0])
        vel = (1.0,)
        m = Microbe{1}(id=1; vel, motility)
        turn!(m, m.motility)
        @test vel == .-m.vel
        turn!(m, m.motility)
        @test vel == m.vel

        motile_state = TwoStates(ForwardState())
        motility = RunReverseFlick(speed_forward=[1.0]; motile_state)
        vel = (1.0,)
        m = Microbe{1}(id=1; vel, motility)
        turn!(m, m.motility)
        @test vel == .-m.vel
        @test motilestate(motility) == BackwardState()
        turn!(m, m.motility)
        @test vel == m.vel
        @test motilestate(motility) == ForwardState()
    end

    @testset "Two-dimensional reorientations" begin
        u = SVector(1.0, 0.0)
        # a null rotation should leave u unchanged
        v1 = rotate(u, 0)
        @test u == v1
        v2 = rotate(u, 0, 0)
        @test u == v2
        # 90-degree rotation
        v3 = rotate(u, π/2)
        @test v3 ≊ SVector(0.0, 1.0)
        # azimuthal angle has no effect on 2d rotation
        v4 = rotate(u, π/2, rand(Uniform(0,2π)))
        @test v3 ≊ v4
        v5 = rotate(u, π)
        # 180-degree rotation
        @test u ≊ .-v5
        # rotation by multiple of 2π leave u unchanged
        v6 = rotate(u, 200π)
        @test u ≊ v6

        U = 30.0
        vel = rand_vel(2) .* U
        motility = RunReverseFlick(
            speed_forward = [U],
            motile_state = TwoStates(ForwardState()),
            polar_backward = [π/2] # exclude -π/2
        )
        m = Microbe{2}(id=1; vel, motility)
        # 5 steps should lead back to initial orientation
        # reverse (π)
        turn!(m, m.motility)
        @test SVector(m.vel) ≊ SVector(-vel[1], -vel[2])
        # flick (π/2)
        turn!(m, m.motility)
        @test dot(m.vel, vel) ≊ 0
        @test SVector(m.vel) ≊ SVector(vel[2], -vel[1])
        # reverse (π)
        turn!(m, m.motility)
        @test SVector(m.vel) ≊ SVector(-vel[2], vel[1])
        # flick (π/2)
        turn!(m, m.motility)
        @test SVector(m.vel) ≊ SVector(-vel[1], -vel[2])
        # reverse (π)
        turn!(m, m.motility)
        @test SVector(m.vel) ≊ SVector(vel)
    end

    @testset "Three-dimensional reorientations" begin
        U = 30.0
        vel = rand_vel(3) .* U
        θ = π/6
        polar = [θ]
        motility = RunTumble(
            speed = [U],
            polar = polar,
        )
        m = Microbe{3}(id=1; vel, motility)
        turn!(m, m.motility)
        @test dot(m.vel, vel)/U^2 ≊ cos(θ)
        
        U = 30.0
        vel = (U, 0.0, 0.0)
        θ = π/4
        φ = π/6
        polar = [θ]
        azimuthal = [φ]
        motility = RunTumble(
            speed = [U],
            polar = polar,
            azimuthal = azimuthal,
        )
        m = Microbe{3}(id=1; vel, motility)
        turn!(m, m.motility)
        @test SVector(m.vel) ≊ SVector(U*cos(θ), U*sin(θ)*sin(φ), -U*sin(θ)*cos(φ))
        @test dot(m.vel, vel)/U^2 ≊ cos(θ)
    end

    @testset "Rotational diffusion" begin
        dt = 1.0
        vel = (1.0,)
        rotational_diffusivity = 0.3
        m = Microbe{1}(id=1; vel, rotational_diffusivity)
        # in 1D rotational diffusion is deactivated, nothing should happen
        rotational_diffusion!(m, dt)
        @test m.vel == vel

        dt = 1.0
        rotational_diffusivity = 0.0
        m = Microbe{2}(id=1; rotational_diffusivity)
        vel = m.vel
        rotational_diffusion!(m, dt)
        # if rotational_diffusivity = 0 vel should be unchanged
        @test m.vel == vel

        dt = 1.0
        rotational_diffusivity = 0.3
        σ = sqrt(2*rotational_diffusivity*dt)
        m = Microbe{2}(id=1; rotational_diffusivity)
        vel = m.vel
        U = norm(vel)
        # fix rng state
        N = 7
        Random.seed!(N)
        θ = rand(Normal(0,σ))
        # reset rng
        Random.seed!(N)
        rotational_diffusion!(m, dt)
        # reorientation should be of an angle θ and norm should be conserved
        @test dot(m.vel, vel)/U^2 ≊ cos(θ)

        # same tests for 3d
        dt = 1.0
        rotational_diffusivity = 0.0
        m = Microbe{3}(id=1; rotational_diffusivity)
        vel = m.vel
        rotational_diffusion!(m, dt)
        # if rotational_diffusivity = 0 vel should be (approx.) unchanged
        @test all(m.vel .≊ vel)

        dt = 1.0
        rotational_diffusivity = 0.3
        σ = sqrt(2*rotational_diffusivity*dt)
        m = Microbe{3}(id=1; rotational_diffusivity)
        vel = m.vel
        U = norm(vel)
        # fix rng state
        N = 7
        Random.seed!(N)
        θ, φ = rand(Normal(0,σ)), rand(Arccos())
        # reset rng
        Random.seed!(N)
        rotational_diffusion!(m, dt)
        # reorientation should be of an angle θ and norm should be conserved
        @test dot(m.vel, vel)/U^2 ≊ cos(θ)
    end
end