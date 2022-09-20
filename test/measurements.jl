using Test, BacteriaBasedModels, Random

@testset "Measurements" begin
    @testset "Chemotactic drift velocity" begin
        m1 = Microbe{1}(id=1, pos=(4.0,), vel=(1.0,), turn_rate=0)
        m2 = Microbe{1}(id=2, pos=(6.0,), vel=(-1.0,), turn_rate=0)
        m3 = Microbe{1}(id=3, pos=(3.0,), vel=(3.0,), turn_rate=0)
        dt = 1.0
        L = 10.0
        target_point = (5.0,)
        target_direction = (1.0,)
        model = initialise_model(;
            microbes = [m1, m2, m3],
            timestep = dt,
            extent = L,
            random_positions = false
        )
        adata = [:pos, :vel]
        adf, = run!(model, microbe_step!, 2; adata)
        adf1 = filter(:id => id -> id==1, adf; view=true)
        adf2 = filter(:id => id -> id==2, adf; view=true)
        adf3 = filter(:id => id -> id==3, adf; view=true)
        vd1_dir = driftvelocity_direction(adf1, target_direction)
        vd1_pnt = driftvelocity_point(adf1, target_point)
        vd2_dir = driftvelocity_direction(adf2, target_direction)
        vd2_pnt = driftvelocity_point(adf2, target_point)
        vd3_pnt = driftvelocity_point(adf3, target_point)
        vd3_pnt_n = driftvelocity_point(adf3, target_point; normalize=true)
        @test vd1_dir == [1.0 1.0 1.0]
        @test vd1_pnt[1] == 1.0
        @test vd1_pnt[3] == -1.0
        @test vd2_dir == [-1.0 -1.0 -1.0]
        @test vd2_pnt[1] == 1.0
        @test vd2_pnt[3] == -1.0
        @test vd3_pnt[1] == 3.0
        @test vd3_pnt[3] == -3.0
        @test vd3_pnt_n[1] == 1.0
        @test vd3_pnt_n[3] == -1.0
    end
end