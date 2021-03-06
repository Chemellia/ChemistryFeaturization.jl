using Documenter, ChemistryFeaturization

makedocs(
	sitename = "ChemistryFeaturization.jl",
	modules = [ChemistryFeaturization],
	pages = Any[
		"Home" => "index.md",
		"Terminology" => "terminology.md",
		"Autodocs" => "autodocs.md"
	]
)
deploydocs(
    repo = "github.com/thazhemadam/ChemistryFeaturization.jl.git",
     target = "build"
)

# deploydocs(
	# repo = "github.com/thazhemadam/ChemistryFeaturization.jl.git",
	# target = "build",
	# push_preview = true,
	# branch = "gh-pages"
# )