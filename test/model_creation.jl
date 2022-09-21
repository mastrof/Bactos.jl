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

    custom_init(extent) = initialise_model(;
        microbes = [Microbe{3}(id=0)],
        extent = extent,
        timestep = dt
    )
    @test_throws ErrorException custom_init((1.0,))
    @test_throws ErrorException custom_init((1.0, 1.0))
    @test_throws ErrorException custom_init((1.0, 1.0, 1.0, 1.0))
    model1 = custom_init(1.0)
    model2 = custom_init((1.0, 2.0, 3.0))
    @test model1.space.extent == (1.0, 1.0, 1.0)
    @test model2.space.extent == (1.0, 2.0, 3.0)
end