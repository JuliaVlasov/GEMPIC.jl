ENV["GKSwstype"] = "100"
using Documenter
using GEMPIC
using Literate
using Plots

makedocs(;
    sitename="GEMPIC",
    doctest=true,
    authors="Julia Vlasov",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", nothing) == "true",
        mathengine=MathJax(
            Dict(
                :TeX =>
                    Dict(:equationNumbers => Dict(:autoNumber => "AMS"), :Macros => Dict()),
            ),
        ),
    ),
    pages=[
        "Documentation" => [
            "index.md",
            "mesh.md",
            "low_level_bsplines.md",
            "splinepp.md",
            "particle_mesh_coupling.md",
            "distributions.md",
            "particle_group.md",
            "particle_sampling.md",
            "maxwell_solver.md",
            "Splitting" => [
                "hamiltonian_splitting.md",
                "hamiltonian_splitting_boris.md",
            ],
            "diagnostics.md",
        ],
        "Examples" => ["strong_landau_damping.md"],
        "Contents" => "contents.md",
    ],
)

deploydocs(; 
    branch = "gh-pages",
    devbranch = "master",
    repo = "github.com/JuliaVlasov/GEMPIC.jl")
