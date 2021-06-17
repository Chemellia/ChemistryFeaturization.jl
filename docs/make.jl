using Documenter, ChemistryFeaturization

makedocs(
    sitename = "ChemistryFeaturization.jl",
    modules = [ChemistryFeaturization],
    pages = Any[
        "Home"=>"index.md",
        "Terminology"=>"terminology.md",
        "Types"=>Any[
            "Atoms Objects"=>"types/atoms.md",
            "Feature Descriptors"=>"types/feature_descriptors.md",
            "Featurization"=>"types/featurizations.md",
        ],
        "Contributing"=>"contributing.md",
        "Changelog"=>"changelog.md",
    ],
    format = Documenter.HTML(
        # Use clean URLs, unless built as a "local" build
        prettyurls = !("local" in ARGS),
        canonical = "https://chemellia.github.io/ChemistryFeaturization.jl/stable/",
        edit_link = "main",
    ),
    linkcheck = "linkcheck" in ARGS,
)
deploydocs(
    repo = "github.com/Chemellia/ChemistryFeaturization.jl.git",
    target = "build",
    devbranch = "main",
    versions = ["stable" => "v^", "v#.#"],
    branch = "gh-pages",
)
