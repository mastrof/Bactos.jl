push!(LOAD_PATH, "../src/")
using Bactos

using Documenter

makedocs(
    sitename = "Bactos.jl",
    modules = [Bactos],
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Validation" => "checks.md"
    ]
)

deploydocs(;
    repo = "github.com/mastrof/Bactos.jl"
)