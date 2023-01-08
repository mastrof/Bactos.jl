using Test, Bactos

@testset "Encounters" begin
    @testset "Spherical targets" begin
        extent = L = 100
        timestep = 0.1
        random_positions = false

        microbes = [Microbe{1}(id=1, turn_rate=0)] # radius=0
        model = initialise_model(; microbes, extent, timestep, random_positions)
        bodies = [ObstacleSphere((0,), 1)]
        add_encounters!(model, bodies)
        # default keys should be added to model
        @test haskey(model.properties, :encounters)
        @test haskey(model.properties, :bodies)
        @test model.encounters == 0
        @test model.bodies === bodies
        # calling add_encounters! again should trigger warnings
        @test_logs (:warn, r"bodies") (:warn, r"encounters") add_encounters!(model, bodies)
        # test keywords
        add_encounters!(model, bodies; key_bodies=:bodies2, key_encounters=:encounters2)
        @test haskey(model.properties, :bodies2)
        @test haskey(model.properties, :encounters2)

        # 2D
        U = 30
        r₁ = 20
        r₂ = 1
        s₁ = ObstacleSphere((L/2,L/2), r₁)
        s₂ = ObstacleSphere((L/2,L/2), r₂)
        function initmodel(s, x₀, v₀=nothing)
            if isnothing(v₀)
                m = Microbe{2}(id=1, pos=x₀, turn_rate=0)
            else
                m = Microbe{2}(id=1, pos=x₀, vel=v₀, turn_rate=0)
            end
            model = initialise_model(;microbes=[m], timestep, extent, random_positions)
            add_encounters!(model, [s])
            model
        end
        v₀ = (U, 0)
        x₀ = (L/2, L/2) .- (v₀.*timestep)
        model1 = initmodel(s₁, x₀, v₀)
        model2 = initmodel(s₂, x₀, v₀)
        # encounters should be triggered with both spheres
        run!(model1) # starts inside s₁ and stays inside s₁
        run!(model2) # starts outside s₂ and moves inside s₂
        @test model1.encounters == 1 && model2.encounters == 1
        # microbes should now have different positions due to random reinsertion
        @test (model1[1].pos ≠ x₀) && (model2[1].pos ≠ x₀) && (model1[1].pos ≠ model2[1].pos)
        # even if the step moves from one side of the sphere to the other,
        # an encounter should be triggered since the displacement intersects the sphere
        s = ObstacleSphere((L/2,L/2), 0.5)
        v₀ = (U, 0) # (30, 0) → during one timestep moves (+3, 0)
        x₀ = @. (L/2, L/2) - v₀*timestep - (0.51, 0) # would cross s with 2nd step
        model = initmodel(s, x₀, v₀)
        run!(model)
        @test model.encounters == 1 && model[1].pos ≠ x₀ .+ (v₀ .* timestep)
    end
end