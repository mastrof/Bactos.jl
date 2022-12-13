using Test, BacteriaBasedModels, Random

@testset "Measurements" begin
    @testset "Chemotactic drift velocity" begin
        # with turn_rate=0 motility is ignored and vel never changes
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
        # m1 is moving along target_direction at all timepoints
        @test vd1_dir == [1.0 1.0 1.0]
        # m1 is moving towards target_point at initial timepoint
        @test vd1_pnt[1] == 1.0
        # at final timepoint, m1 overshoots target_point and moves in opposite direction
        @test vd1_pnt[3] == -1.0
        # m2 is moving against target_direction at all timepoints
        @test vd2_dir == [-1.0 -1.0 -1.0]
        # m2 is moving towards target_point at initial timepoint
        @test vd2_pnt[1] == 1.0
        # at final timepoint, m2 overshoots target_point and moves in opposite direction
        @test vd2_pnt[3] == -1.0
        # same as m1, but now test if normalize kwarg works
        # not normalized
        @test vd3_pnt[1] == 3.0
        @test vd3_pnt[3] == -3.0
        # normalized
        @test vd3_pnt_n[1] == 1.0
        @test vd3_pnt_n[3] == -1.0
    end

    @testset "Mean-squared displacement" begin
        U = 2.0
        extent = 50.0
        # if turn_rate=0 motility is ignored and vel never changes
        # so the microbe is moving ballistically with speed U
        microbes = [Microbe{1}(id=0, turn_rate=0, vel=(U,), pos=(extent/2,))]
        timestep = 0.1
        model = initialise_model(; microbes, extent, timestep,
                                   random_positions=false)
        adata = [:pos]
        nsteps = 10
        adf, = run!(model, microbe_step!, nsteps; adata)
        # total displacement obeys ballistic motion
        Δx = adf[end,:pos][1] .- adf[1,:pos][1]
        @test Δx ≈ U*nsteps*timestep
        MSD = msd(adf) |> vec
        # for ballistic motion MSD(t) = U²t²
        for n in 1:nsteps
            @test MSD[n] ≈ (U*n*timestep)^2
        end # for

        # same ballistic test over longer times
        # this time using `L=extent` keyword to unwrap periodic boundaries
        nsteps = 500
        adf, = run!(model, microbe_step!, nsteps; adata)
        MSD = msd(adf; L=extent) |> vec
        @test MSD[nsteps] ≈ (U*nsteps*timestep)^2

        
        U = 3.0
        # turn_rate=Inf in 1D means that the bacterium is reversing at each step
        microbes = [Microbe{1}(id=0, turn_rate=Inf, vel=(U,),
                    motility=RunTumble(speed=Degenerate(U)))]
        extent = 50.0
        timestep = 0.1 
        model = initialise_model(; microbes, extent, timestep)
        nsteps = 4
        adf, = run!(model, microbe_step!, nsteps; adata)
        MSD = msd(adf; L=extent) |> vec
        # first step is ballistic forward
        @test MSD[1] ≈ (U*timestep)^2
        # second step goes back to initial position
        @test MSD[2] ≈ 0.0
        # same pattern for following timesteps
        @test MSD[3] ≈ MSD[1]
        @test MSD[4] ≈ MSD[2]
    end
end