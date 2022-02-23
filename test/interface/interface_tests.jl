using ChemistryFeaturization.ElementFeature
import ChemistryFeaturization: features

struct DummyFeaturization <: AbstractFeaturization
    features::Vector{ElementFeatureDescriptor}
end
features(df::DummyFeaturization) = df.features

fzn1 =
    DummyFeaturization([ElementFeatureDescriptor("Block"), ElementFeatureDescriptor("X")])

struct DumberFeaturization <: AbstractFeaturization end

fzn2 = DumberFeaturization()

@testset "Interface" begin
    codec_tests = ["feature_tests", "featurization_tests", "featurizedatoms_tests"]
    for t in codec_tests
        tp = abspath(testdir, "interface", "$(t).jl")
        include(tp)
    end
end
