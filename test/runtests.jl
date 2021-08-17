using Conda

# hackaround conda's path woes
ENV["PATH"] = Conda.bin_dir(Conda.ROOTENV) * ";" * ENV["PATH"]

using ChemistryFeaturization
using Test

const testdir = dirname(@__FILE__)

@show ENV["LD_LIBRARY_PATH"]

tests = [
    "module_tests",
    "utils/ElementFeatureUtils_tests",
    "utils/GraphBuilding_tests",
    "atoms/AtomGraph_tests",
    "features/ElementFeature_tests",
    "featurizations/GraphNodeFeaturization_tests",
    "featurizedatoms_tests",
    # TODO: add Weave stuff
    # TODO: add SpeciesFeature tests
]

@testset "ChemistryFeaturization" begin
    for t in tests
        tp = abspath(testdir, "$(t).jl")
        include(tp)
    end
end
