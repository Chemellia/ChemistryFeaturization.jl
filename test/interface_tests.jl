using ChemistryFeaturization.ElementFeature
import ChemistryFeaturization: features

struct DummyFeaturization <: AbstractFeaturization
    features::Vector{ElementFeatureDescriptor}
end
features(df::DummyFeaturization) = df.features

fzn1 = DummyFeaturization([ElementFeatureDescriptor("Block"), ElementFeatureDescriptor("X")])

struct DumberFeaturization <: AbstractFeaturization end

fzn2 = DumberFeaturization()

@testset "Featurization" begin
    @test !("He" in encodable_elements(fzn1))
    @test_throws MethodError features(fzn2)
    
    encoded = encode(C3, fzn1)
    @test size(encoded[1]) == (4,3)
    @test size(encoded[2]) == (10,3)
    @test all(encoded[1][3,:] .== 1)

    decoded = decode(encoded, fzn1)
    @test all(decoded[1] .== "p")
    @test decoded[2][1] == (2.34, 2.668)
end

@testset "FeaturizedAtoms" begin
    fa1 = featurize(C3, fzn1)
    fa2 = FeaturizedAtoms(C3, fzn1)
    @test fa1.encoded_features == fa2.encoded_features
    @test_throws AssertionError featurize(He_mol, fzn1)
    @test decode(fa1) == decode(fa1.encoded_features, fa1.featurization)
end
