using Test, BacteriaBasedModels, Random

@testset "BacteriaBasedModels Tests" begin
    include("microbe_creation.jl")
    include("motility.jl")
    include("model_creation.jl")
    include("reorientations.jl")
    include("step.jl")
    include("measurements.jl")
end