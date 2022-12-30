using Test, Bactos, Random

@testset "Bactos.jl Tests" begin
    include("microbe_creation.jl")
    include("motility.jl")
    include("model_creation.jl")
    include("reorientations.jl")
    include("step.jl")
    include("measurements.jl")
    include("finite_differences.jl")
    include("pathfinder.jl")
    include("encounters.jl")
end