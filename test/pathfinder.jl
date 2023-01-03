using Test, Bactos, Random
using LinearAlgebra: norm

@testset "Pathfinder interface" begin
    walkmap = BitArray(rand((0,1), 10, 10))
    L = 10
    extent = (L, L)
    periodic = false
    pathfinder1 = initialise_pathfinder(extent, periodic, walkmap)
    pathfinder2 = initialise_pathfinder(L, periodic, walkmap)
    # initialising with a tuple or a scalar should be the same
    @test pathfinder1.dims == pathfinder2.dims
    # pathfinder.walkmap should be just a reference to walkmap
    @test pathfinder1.walkmap === walkmap

    spherepos = (5,5)
    sphereradius = 3
    proberadius = 0.5
    spheres = [ObstacleSphere(spherepos, sphereradius)]
    Δ = proberadius # mesh resolution (defaults to proberadius/2)
    walkmap = get_walkmap(extent, proberadius, spheres; Δ)
    # walkmap should be false within the sphere area, true otherwise
    mesh = (x = 0:Δ:L, y = 0:Δ:L)
    walkmap2 = similar(walkmap)
    for j in eachindex(mesh[:y]), i in eachindex(mesh[:x])
        pos = (mesh[:x][i], mesh[:y][j])
        # contact is allowed so if distance = radius the value is true
        if norm(pos .- spherepos) < sphereradius + proberadius
            walkmap2[i,j] = false
        else
            walkmap2[i,j] = true
        end
    end
    @test walkmap == walkmap2

    L = 10
    extent = (L, L)
    timestep = 0.1
    microbes = [Microbe{2}(id=1)]
    model = initialise_model(; microbes, extent, timestep)
    add_pathfinder!(model, walkmap)
    @test haskey(model.properties, :pathfinder)
end