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

    @testset "Differential equations" begin
        microbes = [Microbe{1}(id=0)]
        timestep = 0.1
        extent = 1.0
        model = initialise_model(; microbes, timestep, extent)
        @test !haskey(model.properties, :integrator)
        @test typeof(model.properties) == Dict{Symbol, Float64}
        model = initialise_model(; microbes, timestep, extent, diffeq=true)
        @test haskey(model.properties, :integrator)
        @test model.integrator == BacteriaBasedModels.dummy_integrator
        @test typeof(model.properties) == Dict{Symbol, Any}

        my_ode_step!(du, u, p, t) = (du .= p[1])
        u₀ = [0.0]
        p = (1.0,)
        integrator = initialise_ode(my_ode_step!, u₀, p)
        model = initialise_model(;
            microbes, timestep, extent,
            diffeq = true, ode_integrator = integrator
        )
        @test model.integrator === integrator
        n = 5
        run!(model, microbe_step!, model_step!, n)
        @test model.integrator === integrator
        @test integrator.u[1] ≈ p[1] * timestep * n
    end
end