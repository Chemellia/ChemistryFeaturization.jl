using Documenter, ChemistryFeaturization, PlutoStaticHTML

# this next bit is from https://github.com/rikhuijzer/PlutoStaticHTML.jl/blob/main/docs/make.jl
const NOTEBOOK_DIR = joinpath(pkgdir(ChemistryFeaturization), "docs", "src", "tutorial")

function write_notebooks()
    hopts = HTMLOptions()
    bopts = BuildOptions(NOTEBOOK_DIR; output_format = documenter_output)
    parallel_build(bopts, hopts)
    return nothing
end

write_notebooks()

pages = Any[
    "Overview"=>"index.md",
    "Tutorials"=>Any["Usage: `AtomGraphs`"=>"tutorial/atomgraphs.md",
    # "Implementing"=>"tutorial/interface.md",
    ],
    "API"=>"api.md",
    "Changelog"=>"changelog.md",
]

makedocs(
    sitename = "ChemistryFeaturization.jl",
    modules = [ChemistryFeaturization],
    pages = pages,
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
    branch = "gh-pages",
)

# for local build
"""
makedocs(
    sitename = "ChemistryFeaturization.jl",
    modules = [ChemistryFeaturization],
    pages = pages,
    format = Documenter.HTML(
        # Use clean URLs, unless built as a "local" build
        prettyurls = false,
        canonical = "https://chemellia.github.io/ChemistryFeaturization.jl/stable/",
        edit_link = "main",
    ),
    linkcheck = "linkcheck" in ARGS,
)
"""
