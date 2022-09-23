using Test, BacteriaBasedModels, Random

@testset "Finite differences" begin
    ≃(x,y) = isapprox(x,y;atol=1e-3)
    @testset "One dimension" begin
        n = 256
        x = range(0, 2π; length=n)
        # reset x with ghost cells
        m = 2
        dx = x.step.hi
        x = range(x[1]-m*dx, x[end]+m*dx; step=dx)
        y₀ = sin.(x)
        y₁ = cos.(x)
        y₂ = -sin.(x)
        dy₀ = zero(y₀)
        d²y₀ = zero(y₀)
        ddy₀ = zero(y₀)
        dy₁ = zero(y₀)
        ∇y₀ = zero(y₀)
        ∇²y₀ = zero(y₀)
        divy₀ = zero(y₀)
        finitediff!(dy₀, y₀, 1/dx)
        finitediff!(dy₁, y₁, 1/dx)
        finitediff!(ddy₀, dy₀, 1/dx)
        finitediff!(d²y₀, y₀, 1/dx^2, CFDM_3_2)
        laplacian!(∇y₀, y₀, 1/dx, CFDM_3_1)
        laplacian!(∇²y₀, y₀, 1/dx^2)
        divergence!(divy₀, y₀, 1/dx)
        @test all((dy₀ .≃ y₁)[1+m:end-m])
        @test all((dy₁ .≃ y₂)[1+m:end-m])
        @test all((ddy₀ .≃ y₂)[1+m:end-m])
        @test all((d²y₀ .≃ y₂)[1+m:end-m])
        @test ∇y₀ == dy₀
        @test ∇²y₀ == d²y₀
        @test divy₀ == dy₀
    end

    @testset "Two dimensions" begin
        nx = 256
        ny = 512
        x = range(0, 2π; length=nx)
        y = range(0, 2π; length=ny)
        # reset x and y with ghost cells
        m = 2
        dx = x.step.hi
        dy = y.step.hi
        x = range(x[1]-m*dx, x[end]+m*dx; step=dx)
        y = range(y[1]-m*dy, y[end]-m*dy; step=dy)
        u₀ = @. sin(x) * cos(y)'
        ux = @. cos(x) * cos(y)'
        uy = @. - sin(x) * sin(y)'
        lapu = @. -2 * u₀
        myux = zero(u₀)
        myuy = zero(u₀)
        mylapu = zero(u₀)
        mylapu_2 = zero(u₀)
        finitediff!(myux, myuy, u₀, 1/dx, 1/dy)
        laplacian!(mylapu, u₀, 1/dx^2, 1/dy^2)
        divergence!(mylapu_2, myux, myuy, 1/dx, 1/dy)
        @test all((myux .≃ ux)[1+m:end-m,1+m:end-m])
        @test all((myuy .≃ uy)[1+m:end-m,1+m:end-m])
        @test all((mylapu .≃ lapu)[1+m:end-m,1+m:end-m])
        @test all((mylapu_2 .≃ lapu)[1+m:end-m,1+m:end-m]) # ∇⋅∇ = ∇²
    end
end