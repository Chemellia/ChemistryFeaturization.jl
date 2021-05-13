using ChemistryFeaturization
using Test

const testdir = dirname(@__FILE__)

tests = [
	"utils/AtomFeatureUtils_tests",
	#"atoms/AtomGraph/AtomGraph_tests",
	#"atoms/AtomGraph/graph_building",
	#"atoms/AtomGraph/batch_processing",
	"features/AtomFeature_tests",
	"featurizations/GraphNodeFeaturization_tests",
	# Weave stuff
]

@testset "ChemistryFeaturization" begin
	for t in tests
		tp = abspath(testdir, "$(t).jl")
		include(tp)
	end
end
