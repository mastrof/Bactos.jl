using Test, BacteriaBasedModels, Random
using LinearAlgebra: norm

@testset "Microbe creation" begin
    @testset "Default microbe" begin
        # test type hierarchy
        @test Microbe <: AbstractMicrobe
        @test Microbe{1} <: AbstractMicrobe
        @test Microbe{2} <: AbstractMicrobe
        @test Microbe{3} <: AbstractMicrobe

        # test errors for missing arguments
        @test_throws UndefKeywordError Microbe{3}()
        @test_throws UndefVarError Microbe(id=0)

        # test arguments behave as expected, 2D
        id = rand(Int)
        m = Microbe{2}(id=id)
        @test m.id == id
        @test m.pos isa NTuple{2,Float64}
        @test all(m.pos .== 0.0)
        @test norm(m.vel) ≈ 30.0
        @test typeof(m.motility) == RunTumble
        @test m.turn_rate == 1.0
        @test m.state == 0.0
        @test m.rotational_diffusivity == 0.0
        @test m.radius == 0.0

        # test some arguments again for 3D
        id = rand(Int)
        D = 3
        x = rand(D) .* 500.0 |> Tuple
        v = rand_vel(D) .* 40.0
        mot = RunReverseFlick()
        m = Microbe{D}(id=id, pos=x, vel=v, motility=mot)
        @test m.id == id
        @test m.pos isa NTuple{3,Float64}
        @test m.pos == x
        @test m.vel isa NTuple{3,Float64}
        @test m.vel == v
        @test m.motility isa RunReverseFlick
    end
    
    @testset "Brumley" begin
        # test type hierarchy
        @test MicrobeBrumley{1} <: AbstractMicrobe
        @test MicrobeBrumley{2} <: AbstractMicrobe
        @test MicrobeBrumley{3} <: AbstractMicrobe

        # test arguments
        m = MicrobeBrumley{3}(id=0)
        @test m.pos isa NTuple{3,Float64}
        @test m.vel isa NTuple{3,Float64}
        @test m.turn_rate == 1/0.45
        @test m.state == 0.0
        @test m.rotational_diffusivity == 0.035
        @test m.adaptation_time == 1.3
        @test m.receptor_gain == 50.0
        @test m.motor_gain == 50.0
        @test m.chemotactic_precision == 6.0
        @test m.radius == 0.5
        @test m.motility isa RunReverseFlick
        # default has same speed Degenerate(46.5) for forward and backward modes
        @test norm(m.vel) ≈ rand(m.motility.speed_forward)
        @test norm(m.vel) ≈ rand(m.motility.speed_backward)
    end

    @testset "Brown-Berg" begin
        # test type hierarchy
        @test MicrobeBrownBerg{1} <: AbstractMicrobe
        @test MicrobeBrownBerg{2} <: AbstractMicrobe
        @test MicrobeBrownBerg{3} <: AbstractMicrobe

        # test arguments
        m = MicrobeBrownBerg{3}(id=0)
        @test m.pos isa NTuple{3,Float64}
        @test m.vel isa NTuple{3,Float64}
        @test m.turn_rate == 1/0.67
        @test m.state == 0.0
        @test m.rotational_diffusivity == 0.035
        @test m.motor_gain == 660.0
        @test m.receptor_binding_constant == 100.0
        @test m.adaptation_time == 1.0
        @test m.radius == 0.0
        @test m.motility isa RunTumble
        @test norm(m.vel) ≈ rand(m.motility.speed)
    end

    @testset "Celani" begin
        # test type hierarchy
        @test AbstractCelani{D} <: AbstractMicrobe{D} where {D}
        @test Celani{D} <: AbstractCelani{D} where {D}
        @test CelaniNoisy{D} <: AbstractCelani{D} where {D}
        @test !(CelaniNoisy{D} <: Celani{D} where {D})
        @test !(Celani{D} <: CelaniNoisy{D} where {D})

        @test setdiff(fieldnames(CelaniNoisy), fieldnames(Celani)) == [:chemotactic_precision]

        # test arguments
        m1 = Celani{3}(id=0)
        @test m1.pos isa NTuple{3,Float64}
        @test m1.vel isa NTuple{3,Float64}
        @test m1.turn_rate == 1/0.67
        @test m1.state == [0., 0., 0., 1.]
        @test m1.rotational_diffusivity == 0.26
        @test m1.gain == 50.0
        @test m1.memory == 1.0
        @test m1.radius == 0.5
        @test m1.motility isa RunTumble
        @test norm(m1.vel) ≈ rand(m1.motility.speed)

        m2 = CelaniNoisy{3}(id=0)
        @test m2.pos isa NTuple{3,Float64}
        @test m2.vel isa NTuple{3,Float64}
        @test m2.turn_rate == 1/0.67
        @test m2.state == [0., 0., 0., 1.]
        @test m2.rotational_diffusivity == 0.26
        @test m2.gain == 50.0
        @test m2.memory == 1.0
        @test m2.radius == 0.5
        @test m2.chemotactic_precision == 1.0
        @test m2.motility isa RunTumble
        @test norm(m2.vel) ≈ rand(m2.motility.speed)
    end

    @testset "Xie" begin
        # test type hierarchy
        @test AbstractXie{D} <: AbstractMicrobe{D} where {D}
        @test Xie{D} <: AbstractXie{D} where {D}
        @test XieNoisy{D} <: AbstractXie{D} where {D}
        @test !(XieNoisy{D} <: Xie{D} where {D})
        @test !(Xie{D} <: XieNoisy{D} where {D})

        @test setdiff(fieldnames(XieNoisy), fieldnames(Xie)) == [:chemotactic_precision]

        # test arguments
        m1 = Xie{3}(id=0)
        @test m1.pos isa NTuple{3,Float64}
        @test m1.vel isa NTuple{3,Float64}
        @test m1.turn_rate_forward == 2.3
        @test m1.turn_rate_backward == 1.9
        @test m1.state_m == 0.0
        @test m1.state_z == 0.0
        @test m1.state == 0.0
        @test m1.rotational_diffusivity == 0.26
        @test m1.gain_forward == 2.7
        @test m1.gain_backward == 1.6
        @test m1.adaptation_time_m == 1.29
        @test m1.adaptation_time_z == 0.28
        @test m1.binding_affinity == 0.39
        @test m1.radius == 0.5
        @test m1.motility isa RunReverse
        @test m1.motility.speed_forward == m1.motility.speed_backward

        m2 = XieNoisy{3}(id=0)
        @test m2.pos isa NTuple{3,Float64}
        @test m2.vel isa NTuple{3,Float64}
        @test m2.turn_rate_forward == 2.3
        @test m2.turn_rate_backward == 1.9
        @test m2.state_m == 0.0
        @test m2.state_z == 0.0
        @test m2.state == 0.0
        @test m2.rotational_diffusivity == 0.26
        @test m2.gain_forward == 2.7
        @test m2.gain_backward == 1.6
        @test m2.adaptation_time_m == 1.29
        @test m2.adaptation_time_z == 0.28
        @test m2.binding_affinity == 0.39
        @test m2.radius == 0.5
        @test m2.chemotactic_precision == 6.0
        @test m2.motility isa RunReverse
        @test m2.motility.speed_forward == m2.motility.speed_backward
    end
end