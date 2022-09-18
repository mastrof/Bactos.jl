push!(LOAD_PATH, "../src/")
using BacteriaBasedModels

using Documenter

makedocs(
    sitename = "BacteriaBasedModels.jl",
    modules = [BacteriaBasedModels],
    pages = [
        "Home" => "index.md"
    ]
)

deploydocs(;
    repo = "github.com/mastrof/BacteriaBasedModels"
)