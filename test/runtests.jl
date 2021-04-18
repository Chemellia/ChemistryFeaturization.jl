using ChemistryFeaturization
using Test

const testdir = dirname(@__FILE__)

tests = [
	"AtomFeat/atomfeat",
	"AtomGraph/atomgraph",
	"AtomGraph/batch_processing",
	"AtomGraph/graph_building",
	"graph_visualization/graph_visualization",
	"Weave/smiles",
	"Weave/featurization"
]

@testset "ChemistryFeaturization" begin
	for t in tests
		tp = abspath(testdir, "$(t).jl")
		include(tp)
	end
end
