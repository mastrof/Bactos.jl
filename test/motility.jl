using Test, Bactos, Random
using Distributions: Uniform

@testset "Motility" begin
    @testset "Motile state" begin
        @test instances(TwoState) == (Forward, Backward)
        # call without argument chooses an instance at random
        @test TwoState() ∈ instances(TwoState)
        # bitwise not switches between Forward and Backward
        @test ~Forward == Backward
        @test ~Backward == Forward
        # call to MotileState without argument chooses random state
        motstate = MotileState()
        @test motstate.state ∈ instances(TwoState)
        # or the instance can be specified
        motstate = MotileState(Forward)
        @test motstate.state == Forward
    end

    @test AbstractMotilityOneStep <: AbstractMotility
    @test AbstractMotilityTwoStep <: AbstractMotility
    @test !(AbstractMotilityOneStep <: AbstractMotilityTwoStep)

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
    @test m.motile_state isa MotileState
    # no :state field, but due to getproperty! overload
    # it directly accesses motile_state.state
    @test !hasfield(RunReverse, :state)
    @test m.state === m.motile_state.state
    @test m.state isa TwoState
    # switching motility equals to switching state
    s = m.state
    switch!(m)
    @test m.state == ~s
    # state can be set manually due to setproperty! overload
    m.state = Forward
    @test m.state == Forward

    @test RunReverseFlick <: AbstractMotilityTwoStep
    m = RunReverseFlick()
    @test m.speed_forward == Degenerate(30.0)
    @test m.polar_forward == Degenerate(π)
    @test m.azimuthal_forward == Arccos()
    @test m.speed_backward === m.speed_forward
    @test m.polar_backward == [-π/2, π/2]
    @test m.azimuthal_backward === m.azimuthal_forward
    @test m.motile_state isa MotileState
    # no :state field, but due to getproperty! overload
    # it directly accesses motile_state.state
    @test !hasfield(RunReverse, :state)
    @test m.state === m.motile_state.state
    @test m.state isa TwoState
end