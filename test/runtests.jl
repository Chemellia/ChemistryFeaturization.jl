using ChemistryFeaturization
using Test

const testdir = dirname(@__FILE__)

tests = [
    "interface/interface_tests",
    "ElementFeatureUtils_tests",
    "ElementFeature_tests",
    "codecs/codec_tests",
]

# a few things we'll use in multiple tests
import ChemistryFeaturization: elements
struct DummyAtoms
    els::Vector{String}
end
elements(da::DummyAtoms) = da.els

He_mol = DummyAtoms(["He", "He"])
C3 = DummyAtoms(["C", "C", "C"])

@testset "ChemistryFeaturization" begin
    for t in tests
        tp = abspath(testdir, "$(t).jl")
        include(tp)
    end
end
