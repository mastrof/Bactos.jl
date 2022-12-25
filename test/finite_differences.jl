using Test, Bactos, Random

@testset "Finite differences" begin
    ≃(x,y) = isapprox(x,y;atol=1e-3) # \simeq
    @testset "One dimension" begin
        n = 256
        x = range(0, 2π; length=n)
        # reset x with ghost cells
        m = 2
        dx = x.step.hi
        x = range(x[1]-m*dx, x[end]+m*dx; step=dx)
        # initial function
        y₀ = sin.(x)
        # analytical first derivative
        y₁ = cos.(x)
        # analytical second derivative
        y₂ = -sin.(x)
        # initialize arrays for numerical derivatives
        dy₀ = zero(y₀)
        d²y₀ = zero(y₀)
        ddy₀ = zero(y₀)
        dy₁ = zero(y₀)
        ∇y₀ = zero(y₀)
        ∇²y₀ = zero(y₀)
        divy₀ = zero(y₀)
        # evaluate numerical derivatives
        # first derivative
        finitediff!(dy₀, y₀, 1/dx)
        # first derivative of analytical first derivative
        finitediff!(dy₁, y₁, 1/dx)
        # first derivative of numerical first derivative
        finitediff!(ddy₀, dy₀, 1/dx)
        # second derivative
        finitediff!(d²y₀, y₀, 1/dx^2, CFDM_3_2)
        # with a 1st order stencil, laplacian equates first derivative
        laplacian!(∇y₀, y₀, 1/dx, CFDM_3_1)
        # actual laplacian, in 1d corresponds to second derivative
        laplacian!(∇²y₀, y₀, 1/dx^2)
        # divergence, in 1d corresponds to first derivative
        divergence!(divy₀, y₀, 1/dx)
        # ghost cells must be excluded from the test
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

    @testset "Three dimensions" begin
        nx = 128
        ny = 128
        nz = 64
        x₁, x₂ = 0.0, 2π
        y₁, y₂ = 0.0, 2π
        z₁, z₂ = 0.0, 1.0
        dx = (x₂-x₁)/(nx-1)
        dy = (y₂-y₁)/(ny-1)
        dz = (z₂-z₁)/(nz-1)
        # domain with ghost cells
        m = 2
        x = range(x₁-m*dx, x₂+m*dx; step=dx)
        y = range(y₁-m*dy, y₂+m*dy; step=dy)
        z = range(z₁-m*dz, z₂+m*dz; step=dz)
        u₀ = zeros(nx+2m, ny+2m, nz+2m)
        ux = zero(u₀)
        uy = zero(u₀)
        uz = zero(u₀)
        ∇²u = zero(u₀)
        for k in axes(u₀,3), j in axes(u₀,2), i in axes(u₀,1)
            u₀[i,j,k] = cos(x[i])*sin(y[j])*z[k]
            ux[i,j,k] = -sin(x[i])*sin(y[j])*z[k]
            uy[i,j,k] = cos(x[i])*cos(y[j])*z[k]
            uz[i,j,k] = cos(x[i])*sin(y[j])
            ∇²u[i,j,k] = -2*u₀[i,j,k]
        end
        myux = zero(u₀)
        myuy = zero(u₀)
        myuz = zero(u₀)
        my∇²u = zero(u₀)
        finitediff!(myux, myuy, myuz, u₀, 1/dx, 1/dy, 1/dz)
        laplacian!(my∇²u, u₀, 1/dx^2, 1/dy^2, 1/dz^2)
        @test all((myux .≃ ux)[m+1:end-m,m+1:end-m,m+1:end-m])
        @test all((myuy .≃ uy)[m+1:end-m,m+1:end-m,m+1:end-m])
        @test all((myuz .≃ uz)[m+1:end-m,m+1:end-m,m+1:end-m])
        @test all((my∇²u .≃ ∇²u)[m+1:end-m,m+1:end-m,m+1:end-m])
    end
end