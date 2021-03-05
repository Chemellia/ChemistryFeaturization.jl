using Documenter, ChemistryFeaturization

makedocs(
	sitename = "ChemistryFeaturization.jl",
	modules = [ChemistryFeaturization],
	pages = Any[
		"Home" => "index.md",
	]
)

# deploydocs(
	# repo = "github.com/thazhemadam/ChemistryFeaturization.jl.git",
	# target = "build",
	# push_preview = true,
	# branch = "gh-pages"
# )