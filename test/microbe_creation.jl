using Test, BacteriaBasedModels, Random
using LinearAlgebra: norm

@testset "Microbe creation" begin
    @testset "Default microbe" begin
        @test Microbe <: AbstractMicrobe
        @test Microbe{1} <: AbstractMicrobe
        @test Microbe{2} <: AbstractMicrobe
        @test Microbe{3} <: AbstractMicrobe

        id = rand(Int)
        m = Microbe{2}(id=id)
        @test m.id == id
        @test all(m.pos .== 0.0)
        @test norm(m.vel) ≈ 1.0
        @test typeof(m.motility) == RunTumble

        id = rand(Int)
        D = 3
        x = rand(D) .* 500.0 |> Tuple
        v = rand_vel(D) .* 40.0
        mot = RunReverseFlick()
        m = Microbe{D}(id=id, pos=x, vel=v, motility=mot)
        @test m.id == id
        @test m.pos == x
        @test m.vel == v
        @test typeof(m.motility) == RunReverseFlick
    end
    
    @testset "Brumley" begin
        @test MicrobeBrumley{1} <: AbstractMicrobe
        @test MicrobeBrumley{2} <: AbstractMicrobe
        @test MicrobeBrumley{3} <: AbstractMicrobe

        m = MicrobeBrumley{3}(id=0)
        @test typeof(m.motility) == RunReverseFlick
        @test sqrt(sum(abs2.(m.vel))) ≈ rand(m.motility.speed)

        v = 1.0
        mot = RunReverse(speed=Degenerate(v))
        m = MicrobeBrumley{3}(id=0, motility=mot)
        @test typeof(m.motility) == RunReverse
        @test sqrt(sum(abs2.(m.vel))) ≈ v
    end

    @testset "Brown-Berg" begin
        @test MicrobeBrownBerg{1} <: AbstractMicrobe
        @test MicrobeBrownBerg{2} <: AbstractMicrobe
        @test MicrobeBrownBerg{3} <: AbstractMicrobe

        m = MicrobeBrownBerg{3}(id=0)
        @test typeof(m.motility) == RunTumble
        @test sqrt(sum(abs2.(m.vel))) ≈ rand(m.motility.speed)

        v = 1.0
        mot = RunReverse(speed=Degenerate(v))
        m = MicrobeBrownBerg{3}(id=0, motility=mot)
        @test typeof(m.motility) == RunReverse
        @test sqrt(sum(abs2.(m.vel))) ≈ v
    end
end