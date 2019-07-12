using Documenter
using GEMPIC

makedocs(
    sitename = "GEMPIC",
    format = Documenter.HTML(),
    modules = [GEMPIC],
    pages = ["Documentation" => ["index.md",
                                 "mesh.md",
                                 "low_level_bsplines.md",
                                 "splinepp.md",
                                 "particle_mesh_coupling.md",
                                 "distributions.md",
                                 "particle_group.md",
                                 "particle_sampling.md",
                                 "maxwell_solver.md",
                                 "hamiltonian_splitting.md",
                                 "hamiltonian_splitting_boris.md",
                                 "diagnostics.md"],
             "Contents"      => "contents.md"]
)

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/JuliaVlasov/GEMPIC.jl.git",
 )
