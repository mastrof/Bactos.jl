using Test, Bactos, Random

@testset "Model creation" begin
    timestep = 1.0
    extent = 500.0
    pos = Tuple(rand(3))
    m = Microbe{3}(id=1, pos = pos)
    microbes = [m]

    # test errors
    # `microbes`, `timestep` and `extent` MUST all be given
    @test_throws UndefKeywordError initialise_model()
    @test_throws UndefKeywordError initialise_model(;microbes)
    @test_throws UndefKeywordError initialise_model(;timestep)
    @test_throws UndefKeywordError initialise_model(;extent)
    @test_throws UndefKeywordError initialise_model(;microbes,timestep)
    @test_throws UndefKeywordError initialise_model(;microbes,extent)
    @test_throws UndefKeywordError initialise_model(;timestep,extent)

    model = initialise_model(;
        microbes,
        timestep,
        extent,
        random_positions = false,
    )
    # test default properties
    @test model.properties isa Dict{Symbol,Any}
    @test Set(keys(model.properties)) == Set(
        (:timestep, :compound_diffusivity, :concentration_field,
        :concentration_gradient, :concentration_time_derivative, :integrator)
    )
    @test model.timestep == timestep
    # the model should contain the agent `m`, not a copy
    @test model.agents[1] === m
    # agent position should be conserved since random_positions=false
    @test model.agents[1].pos == pos

    model = initialise_model(;
        microbes,
        timestep,
        extent,
        random_positions = true,
    )
    # agent identity is still the same
    @test model.agents[1] === m
    # agent position is changed since random_positions=true
    @test model.agents[1].pos ≠ pos
    @test m.pos ≠ pos
    # the new position should be inside of extent
    @test all(0 .≤ model.agents[1].pos .< extent)

    # if extent is a tuple, orthorhombic domains can be created
    model = initialise_model(;
        microbes = [Microbe{3}(id=0)], timestep,
        extent = (300.0,400.0,250.0)
    )
    @test model.space.extent == (300.0, 400.0, 250.0)

    # if extent is a scalar, the domain is cubic;
    # if extent is a tuple, the domain can be orthorhombic;
    # if extent has different size from microbe dimensionality,
    # errors should be thrown
    my_init(extent) = initialise_model(;
        microbes = [Microbe{3}(id=0)],
        extent,
        timestep
    )
    @test_throws ErrorException my_init((1.0,))
    @test_throws ErrorException my_init((1.0, 1.0))
    @test_throws ErrorException my_init((1.0, 1.0, 1.0, 1.0))
    model1 = my_init(1.0)
    model2 = my_init((1.0, 2.0, 3.0))
    @test model1.space.extent == (1.0, 1.0, 1.0)
    @test model2.space.extent == (1.0, 2.0, 3.0)

    # test mixed species models
    timestep = 1.0
    extent = 1.0
    microbes = [MicrobeBrumley{1}(id=1), MicrobeBrumley{1}(id=2)]
    model = initialise_model(; microbes, timestep, extent)
    @test model.agents isa Dict{Int, MicrobeBrumley{1}}
    microbes = [Microbe{1}(id=1), MicrobeBrumley{1}(id=2), Xie{1}(id=3)]
    model = initialise_model(; microbes, timestep, extent)
    @test model.agents isa Dict{Int, Union{Microbe{1}, MicrobeBrumley{1}, Xie{1}}}

    @testset "DiffEq Integrator" begin
        microbes = [Microbe{1}(id=0)]
        timestep = 0.1
        extent = 1.0
        model = initialise_model(; microbes, timestep, extent)

        # simple ode: du/dt = p
        # solution: u(t) = u(0) + p*t
        my_ode_step!(du, u, p, t) = (du .= p[1])
        u₀ = [0.0]
        p = (1.0,)
        integrator = initialise_ode(my_ode_step!, u₀, p)
        # check integrator is initialised correctly
        @test integrator.u == u₀
        @test !(integrator.u === u₀)
        @test integrator.p == p
        @test integrator.p === p

        model = initialise_model(;
            microbes, timestep, extent,
            ode_integrator = integrator
        )
        @test model.integrator === integrator
        # advance ode for n steps of size timestep
        n = 5
        run!(model, microbe_step!, model_step!, n)
        @test model.integrator === integrator
        @test integrator.u[1] ≈ p[1] * timestep * n
    end
end