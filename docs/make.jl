push!(LOAD_PATH, "../src/")
ENV["GKSwstype"]="100"
using Documenter
using GEMPIC
using Literate
using Plots

# generate examples
output = joinpath(@__DIR__, "src", "generated")
examples = String[]
push!(examples, "strong_landau_damping.jl")

for example in examples
    jl_file = joinpath(@__DIR__, "examples", example)
    Literate.markdown(jl_file, output, documenter=true)
    #Literate.notebook(EXAMPLE, OUTPUT)
    #Literate.script(EXAMPLE, OUTPUT)
end

makedocs(
    sitename = "GEMPIC",
    doctest = true,
    authors = "Julia Vlasov", 
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        mathengine = MathJax(Dict(
            :TeX => Dict(
                :equationNumbers => Dict(:autoNumber => "AMS"),
                :Macros => Dict()
            )
        ))
    ),
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
                                 "Splitting" => [
                                 "hamiltonian_splitting.md",
                                 "hamiltonian_splitting_spin.md",
                                 "hamiltonian_splitting_boris.md"],
                                 "diagnostics.md"],
              "Examples" => [ "generated/strong_landau_damping.md" ],
             # "Scalar Spin Vlasov-Maxwell" => "scalar_spin_vlasov_maxwell.md",
             "Contents"      => "contents.md"]
)

deploydocs(
    repo   = "github.com/JuliaVlasov/GEMPIC.jl.git",
	)
