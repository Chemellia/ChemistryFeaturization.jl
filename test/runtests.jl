using ChemistryFeaturization
using Test

const testdir = dirname(@__FILE__)

tests = [
	"utils/AtomFeatureUtils_tests",
	"utils/GraphBuilding_tests",
	"atoms/AtomGraph_tests",
	"features/AtomFeature_tests",
	"featurizations/GraphNodeFeaturization_tests",
	# TODO: add Weave stuff
]

@testset "ChemistryFeaturization" begin
	for t in tests
		tp = abspath(testdir, "$(t).jl")
		include(tp)
	end
end
