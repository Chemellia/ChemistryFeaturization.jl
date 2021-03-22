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
    ],
	format = Documenter.HTML(
        # Use clean URLs, unless built as a "local" build
        prettyurls = !("local" in ARGS),
        canonical = "https://chemellia.github.io/ChemistryFeaturization.jl/stable/",
        edit_link = "main",
    ),
    linkcheck = "linkcheck" in ARGS
)
deploydocs(
	repo = "github.com/thazhemadam/ChemistryFeaturization.jl.git",
        target = "build",
        devbranch = "main",
	branch = "gh-pages"
)
