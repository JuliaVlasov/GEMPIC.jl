using Documenter
using GEMPIC

makedocs(
    sitename = "GEMPIC",
    format = Documenter.HTML(),
    modules = [GEMPIC]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
