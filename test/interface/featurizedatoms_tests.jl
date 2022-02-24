@testset "FeaturizedAtoms" begin
    fa1 = featurize(C3, fzn1)
    fa2 = FeaturizedAtoms(C3, fzn1)
    @test fa1.encoded_features == fa2.encoded_features
    @test_throws AssertionError featurize(He_mol, fzn1)
    @test decode(fa1) == decode(fa1.encoded_features, fa1.featurization)
end
