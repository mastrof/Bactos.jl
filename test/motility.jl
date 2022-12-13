using Test, BacteriaBasedModels, Random
using Distributions: Uniform

@testset "Motility" begin
    @test AbstractMotilityOneStep <: AbstractMotility
    @test AbstractMotilityTwoStep <: AbstractMotility
    @test !(AbstractMotilityOneStep <: AbstractMotilityTwoStep)
    @test TwoStates <: AbstractMotileState

    @test RunTumble <: AbstractMotilityOneStep
    m = RunTumble()
    @test m.speed == Degenerate(30.0)
    @test m.polar == Uniform(-π, π)
    @test m.azimuthal == Arccos()

    @test RunReverse <: AbstractMotilityTwoStep
    m = RunReverse()
    @test m.speed_forward == Degenerate(30.0)
    @test m.polar_forward == Degenerate(π)
    @test m.azimuthal_forward == Arccos()
    @test m.speed_backward === m.speed_forward
    @test m.polar_backward === m.polar_forward
    @test m.azimuthal_backward === m.azimuthal_forward
    @test m.motile_state isa TwoStates
    @test motilestate(m) ∈ (ForwardState(), BackwardState())

    @test RunReverseFlick <: AbstractMotilityTwoStep
    m = RunReverseFlick()
    @test m.speed_forward == Degenerate(30.0)
    @test m.polar_forward == Degenerate(π)
    @test m.azimuthal_forward == Arccos()
    @test m.speed_backward === m.speed_forward
    @test m.polar_backward == Degenerate(π/2)
    @test m.azimuthal_backward === m.azimuthal_forward
    @test m.motile_state isa TwoStates
    @test motilestate(m) ∈ (ForwardState(), BackwardState())
end