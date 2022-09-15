using Test, BacteriaBasedModels, Random

@testset "Model creation" begin
    dt = 0.1
    pos = Tuple(rand(3))
    m = Microbe{3}(id=1, pos = pos)
    model = initialise_model(;
        microbes = [m],
        extent = 1.0,
        random_positions = false,
        timestep = dt
    )
    @test model.properties == Dict(:timestep => dt)
    @test model.timestep == dt
    @test model.agents[1] === m
    @test model.agents[1].pos == pos

    model = initialise_model(;
        microbes = [m],
        extent = 1.0,
        random_positions = true,
        timestep = dt
    )
    @test model.agents[1] === m
    @test model.agents[1].pos ≠ pos
    @test all(model.agents[1].pos .< 1)
end