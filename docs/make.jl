using Documenter
using GEMPIC

makedocs(
    sitename = "GEMPIC",
    format = Documenter.HTML(),
    modules = [GEMPIC]
)

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/juliavlasov/GEMPIC.jl.git",
 )
