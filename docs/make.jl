using Documenter, ChemistryFeaturization

makedocs(
	sitename = "ChemistryFeaturization.jl",
	modules = [ChemistryFeaturization],
	pages = Any[
		"Home" => "index.md",
		"Terminology" => "terminology.md",
	]
)
deploydocs(
	repo = "github.com/thazhemadam/ChemistryFeaturization.jl.git",
	target = "build",
	branch = "gh-pages"
)
