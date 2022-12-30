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
        microbes = [Microbe{2}(id=1, pos=(L/2,L/2,L/2))] # radius=0
        model = initialise_model(; microbes, extent, timestep, random_positions)
        @test is_encounter(microbes[1], s₁, model)
        @test is_encounter(microbes[1], s₂, model)
        # if microbe starts from surface contact, it always triggers an encounter
        microbes = [Microbe{2}(id=1, pos=(L/2-r₁,L/2,L/2))]
        model = initialise_model(; microbes, extent, timestep, random_positions)
        @test is_encounter(microbes[1], s₁, model)
        # if a microbe reaches contact distance, it triggers an encounter
        microbes = [Microbe{2}(id=1, pos=(L/2-r₁-U*timestep,L/2,L/2), vel=(U,0,0))]
        model = initialise_model(; microbes, extent, timestep, random_positions)
        @test is_encounter(microbes[1], s₁, model)
        # encounter triggered even if *both* endpoints are outside sphere
        microbes = [Microbe{2}(id=1, pos=(L/2-r₂-0.5,L/2,L/2), vel=(U,0,0))]
        model = initialise_model(; microbes, extent, timestep, random_positions)
        @test is_encounter(microbes[1], s₂, model)
        # if displacement does not cross the sphere, no encounter is triggered
        microbes = [Microbe{2}(id=1, pos=(L/2-2r₁,L/2,L/2), vel=(0,U,0))]
        model = initialise_model(; microbes, extent, timestep, random_positions)
        @test ~is_encounter(microbes[1], s₁, model)
        @test ~is_encounter(microbes[1], s₁, model)
    end
end