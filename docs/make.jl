using Documenter

push!(LOAD_PATH, "../src/")
using TwoDimensional

DEPLOYDOCS = (get(ENV, "CI", nothing) == "true")

makedocs(
    sitename = "TwoDimensional.jl Package",
    format = Documenter.HTML(
        prettyurls = DEPLOYDOCS,
    ),
    authors = "Éric Thiébaut and contributors",
    pages = ["index.md", "install.md", "AffineTransforms.md",
             "Points.md", "BoundingBoxes.md", "reference.md"]
)

if DEPLOYDOCS
    deploydocs(
        repo = "github.com/emmt/TwoDimensional.jl.git",
    )
end
