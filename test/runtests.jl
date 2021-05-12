using ChemistryFeaturization
using Test

const testdir = dirname(@__FILE__)

tests = [
	"AtomGraph/AtomGraph_tests",
	"AtomGraph/batch_processing",
	"AtomGraph/graph_building",
	"AtomFeature/AtomFeature_tests",
	"AtomFeature/AtomFeatureUtils_tests",
	# featurization
	"Weave/smiles",
	"Weave/featurization"
]

@testset "ChemistryFeaturization" begin
	for t in tests
		tp = abspath(testdir, "$(t).jl")
		include(tp)
	end
end
