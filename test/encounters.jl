using Test, Bactos

@testset "Encounters" begin
    @testset "Spherical targets" begin
        extent = L = 100
        timestep = 0.1
        random_positions = false

        # 2D
        x₀ = (L/2, L/2)
        U = 30
        r₁ = 20
        r₂ = 1
        s₁ = ObstacleSphere(x₀, r₁)
        s₂ = ObstacleSphere(x₀, r₂)
        # microbe inside the sphere always triggers an encounter
        microbes = [Microbe{2}(id=1, pos=(L/2,L/2))] # radius=0
        model = initialise_model(; microbes, extent, timestep, random_positions)
        @test is_encounter(microbes[1], s₁, model)
        @test is_encounter(microbes[1], s₂, model)
        # if microbe starts from surface contact, it always triggers an encounter
        microbes = [Microbe{2}(id=1, pos=(L/2-r₁,L/2))]
        model = initialise_model(; microbes, extent, timestep, random_positions)
        @test is_encounter(microbes[1], s₁, model)
        # if a microbe reaches contact distance, it triggers an encounter
        microbes = [Microbe{2}(id=1, pos=(L/2-r₁-U*timestep,L/2), vel=(U,0))]
        model = initialise_model(; microbes, extent, timestep, random_positions)
        @test is_encounter(microbes[1], s₁, model)
        # encounter triggered even if *both* endpoints are outside sphere
        microbes = [Microbe{2}(id=1, pos=(L/2-r₂-0.5,L/2), vel=(U,0))]
        model = initialise_model(; microbes, extent, timestep, random_positions)
        @test is_encounter(microbes[1], s₂, model)
        # if displacement does not cross the sphere, no encounter is triggered
        microbes = [Microbe{2}(id=1, pos=(L/2-2r₁,L/2), vel=(0,U))]
        model = initialise_model(; microbes, extent, timestep, random_positions)
        @test ~is_encounter(microbes[1], s₁, model)
        @test ~is_encounter(microbes[1], s₁, model)

        # repeat all tests in 3D just to be sure
        x₀ = (L/2, L/2, L/2)
        U = 30
        r₁ = 20
        r₂ = 1
        s₁ = ObstacleSphere(x₀, r₁)
        s₂ = ObstacleSphere(x₀, r₂)
        # microbe inside the sphere always triggers an encounter
        microbes = [Microbe{3}(id=1, pos=(L/2,L/2,L/2))] # radius=0
        model = initialise_model(; microbes, extent, timestep, random_positions)
        @test is_encounter(microbes[1], s₁, model)
        @test is_encounter(microbes[1], s₂, model)
        # if microbe starts from surface contact, it always triggers an encounter
        microbes = [Microbe{3}(id=1, pos=(L/2-r₁,L/2,L/2))]
        model = initialise_model(; microbes, extent, timestep, random_positions)
        @test is_encounter(microbes[1], s₁, model)
        # if a microbe reaches contact distance, it triggers an encounter
        microbes = [Microbe{3}(id=1, pos=(L/2-r₁-U*timestep,L/2,L/2), vel=(U,0,0))]
        model = initialise_model(; microbes, extent, timestep, random_positions)
        @test is_encounter(microbes[1], s₁, model)
        # encounter triggered even if *both* endpoints are outside sphere
        microbes = [Microbe{3}(id=1, pos=(L/2-r₂-0.5,L/2,L/2), vel=(U,0,0))]
        model = initialise_model(; microbes, extent, timestep, random_positions)
        @test is_encounter(microbes[1], s₂, model)
        # if displacement does not cross the sphere, no encounter is triggered
        microbes = [Microbe{3}(id=1, pos=(L/2-2r₁,L/2,L/2), vel=(0,U,0))]
        model = initialise_model(; microbes, extent, timestep, random_positions)
        @test ~is_encounter(microbes[1], s₁, model)
        @test ~is_encounter(microbes[1], s₁, model)

        # test the ABM interface
        extent = L = 100
        timestep = 0.1
        random_positions = false
        U = 30
        r = 20
        s = ObstacleSphere((L/2,L/2), r)
        spheres = [s]
        microbes = [Microbe{2}(id=1, pos=(L/2-r-U*timestep/2,L/2), vel=(U,0))]
        x₀ = microbes[1].pos
        model_properties = Dict(:bodies => spheres)
        model = initialise_model(;
            microbes, timestep, extent, random_positions,
            model_properties
        )
        my_model_step!(model) = model_step!(model, update_model! = encounters!)
        # an encounter should occur at first step
        run!(model, microbe_step!, my_model_step!, 1)
        @test model.encounters == 1
        # and the microbe should be in a new position
        @test ~all(model[1].pos .≈ x₀ .+ (U,0).*timestep)
        # the new position is also not triggering an encounter in the next step
        run!(model, microbe_step!, my_model_step!, 1)
        @test model.encounters == 1
        # test reinsert! explicitly
        x₁ = model[1].pos
        reinsert!(model[1], model, :bodies)
        # new position which does not trigger an encounter
        @test ~all(model[1].pos .≈ x₁)
        run!(model, microbe_step!, my_model_step!, 1)
        @test model.encounters == 1
    end
end