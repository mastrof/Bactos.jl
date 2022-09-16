using Test, BacteriaBasedModels, Random
using LinearAlgebra: norm
using Agents: step!

@testset "Stepping" begin
    @testset "Microbe stepping" begin
        pos = ntuple(_ -> 0.5, 2)
        vel = (-0.2, 0.6)
        m = Microbe{2}(id=1, pos=pos, vel=vel, turn_rate=0.0, rotational_diffusivity=0.0)
        dt = 0.1
        model = initialise_model(;
            timestep = dt,
            microbes = [m],
            extent = 1.0, periodic = true,
            random_positions = false
        )

        step!(model, microbe_step!)
        @test m.pos[1] ≈ 0.48
        @test m.pos[2] ≈ 0.56
        @test m.vel == vel

        step!(model, microbe_step!, 9)
        @test m.pos[1] ≈ 0.3
        @test m.pos[2] ≈ 0.1 # periodic boundary conditions!
        @test m.vel == vel
    end
end