using ChemistryFeaturization
using Test

const testdir = dirname(@__FILE__)

tests = [
    "module_tests",
    "utils/ElementFeatureUtils_tests",
    "utils/GraphBuilding_tests",
    "atoms/AtomGraph_tests",
    "features/ElementFeature_tests",
    "features/OrbitalFeature_tests",
    "utils/OrbitalFeatureUtils_tests",
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
