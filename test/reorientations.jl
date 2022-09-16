using Test, BacteriaBasedModels, Random
using LinearAlgebra: norm, dot
using StaticArrays
using Distributions

@testset "Microbe reorientations" begin
    ≊(x,y) = isapprox(x, y, atol=1e-12)

    @testset "One-dimensional reorientations" begin
        vel = (1.0,)
        m = Microbe{1}(id=1, vel=vel, motility=RunTumble())
        turn!(m, m.motility)
        @test vel == .-m.vel
        turn!(m, m.motility)
        @test vel == m.vel

        m = Microbe{1}(id=1, vel=vel, motility=RunReverseFlick())
        s = m.motility.motile_state[1]
        turn!(m, m.motility)
        @test vel == .-m.vel
        @test s == 1 - m.motility.motile_state[1]
        turn!(m, m.motility)
        @test vel == m.vel
        @test s == m.motility.motile_state[1]
    end

    @testset "Two-dimensional reorientations" begin
        u = SVector(1.0, 0.0)
        v1 = rotate(u, 0)
        @test u == v1
        v2 = rotate(u, 0, 0)
        @test u == v2
        v3 = rotate(u, π/2)
        @test dot(u, v3) ≊ 0
        v4 = rotate(u, π/2, rand(Uniform(0,2π)))
        @test dot(u, v4) ≊ 0
        @test all(v3 .≊ v4)
        v5 = rotate(u, π)
        @test u ≊ .-v5
        v6 = rotate(u, 200π)
        @test u ≊ v6

        vel = rand_vel(2)
        motility = RunReverseFlick(
            yaw_0 = Degenerate(π),
            yaw_1 = Degenerate(π/2),
            motile_state = [0]
        )
        m = Microbe{2}(id=1, vel=vel, motility=motility)
        turn!(m, m.motility)
        @test dot(m.vel, vel) ≊ -1
        @test m.motility.motile_state[1] == 1
        turn!(m, m.motility)
        @test dot(m.vel, vel) ≊ 0
        @test m.motility.motile_state[1] == 0
        turn!(m, m.motility)
        @test dot(m.vel, vel) ≊ 0
        turn!(m, m.motility)
        @test dot(m.vel, vel) ≊ -1
        turn!(m, m.motility)
        @test dot(m.vel, vel) ≊ 1
    end

    @testset "Three-dimensional reorientations" begin
        vel = rand_vel(3)
        motility = RunReverseFlick(
            yaw_0 = Degenerate(π),
            pitch_0 = Uniform(0,2π),
            yaw_1 = Degenerate(π/2),
            pitch_1 = Uniform(0,2π),
            motile_state = [0]
        )
        m = Microbe{3}(id=1, vel=vel, motility=motility)
        turn!(m, m.motility)
        @test dot(m.vel, vel) ≊ -1
        turn!(m, m.motility)
        @test dot(m.vel, vel) ≊ 0
        turn!(m, m.motility)
        @test dot(m.vel, vel) ≊ 0
    end
end