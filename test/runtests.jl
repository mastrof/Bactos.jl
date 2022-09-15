using Test, BacteriaBasedModels, Random

@testset "BacteriaBasedModels Tests" begin
    include("microbe_creation.jl")
    include("model_creation.jl")
end