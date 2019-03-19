using GEMPIC
using Documenter

makedocs(modules=[GEMPIC],
         doctest = false,
         format = :html,
         sitename = "GEMPIC.jl",
         pages = ["Home" => "index.md"],
                  "Contents" => "contents.md"])

deploydocs(
    deps   = Deps.pip("mkdocs", "python-markdown-math"),
    repo   = "github.com:JuliaVlasov/GEMPIC.jl.git",
 )
