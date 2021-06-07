using Documenter, ChemistryFeaturization

makedocs(
    sitename = "ChemistryFeaturization.jl",
    modules = [ChemistryFeaturization],
    pages = Any[
    "Home" => "index.md",
    "Terminology" => "terminology.md",
    "AtomFeat" => "types/atomfeat.md",
    "AtomGraph" => "types/atomgraph.md",
    "Building Atomic Graphs" => "build_ag.md",
        "Changelog" => "changelog.md",    
    ]
)
deploydocs(
    repo = "github.com/Chemellia/ChemistryFeaturization.jl.git",
    target = "build",
    devbranch = "main",
    branch = "gh-pages"
)
