using Documenter, ChemistryFeaturization

makedocs(
    sitename = "ChemistryFeaturization.jl",
    modules = [ChemistryFeaturization],
    pages = Any[
        "Home"=>"index.md",
        "Terminology"=>"terminology.md",
        "Tutorial"=>"tutorial.md",
        "Types"=>Any[
            "Overview" => "types/overview.md",
            "Abstract Types"=>"types/abstracttypes.md",
            "Atoms Objects"=>"types/atoms.md",
            "Feature Descriptors"=>"types/feature_descriptors.md",
            "Codec"=>"types/codecs.md",
            "Featurization"=>"types/featurizations.md",
            "Featurized Atoms" => "types/featurizedatoms.md"
        ],
        "Utilities"=>"utils.md",
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

# for local build
"""
makedocs(
    sitename = "ChemistryFeaturization.jl",
    modules = [ChemistryFeaturization],
    pages = Any[
        "Home"=>"index.md",
        "Terminology"=>"terminology.md",
        "Tutorial"=>"tutorial.md",
        "Types"=>Any[
            "Abstract Types"=>"types/abstracttypes.md",
            "Atoms Objects"=>"types/atoms.md",
            "Feature Descriptors"=>"types/feature_descriptors.md",
            "Featurization"=>"types/featurizations.md",
            "Featurized Atoms"=>"types/featurizedatoms.md",
            "Codec"=>"types/codecs.md",
        ],
        "Utilities"=>"utils.md",
        "Contributing"=>"contributing.md",
        "Changelog"=>"changelog.md",
    ],
    format = Documenter.HTML(
        # Use clean URLs, unless built as a "local" build
        prettyurls = false,
        canonical = "https://chemellia.github.io/ChemistryFeaturization.jl/stable/",
        edit_link = "main",
    ),
    linkcheck = "linkcheck" in ARGS,
)
"""
