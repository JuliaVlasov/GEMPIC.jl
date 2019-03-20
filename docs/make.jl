using Documenter
using GEMPIC

makedocs(
    sitename = "GEMPIC",
    format = Documenter.HTML(),
    modules = [GEMPIC]
    pages = ["Documentation" => "index.md",
             "Contents"      => "contents.md"]
)

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com/JuliaVlasov/GEMPIC.jl.git",
 )
