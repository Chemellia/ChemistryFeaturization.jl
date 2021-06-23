using Test
using ChemistryFeaturization.Featurization
using ChemistryFeaturization.Atoms: AtomGraph

@testset "DummyFeaturization" begin
    df = DummyFeaturization()
    @test isempty(encodable_elements(df))
    dummy_encoded = [0 1 0]
    @test_throws ArgumentError decode(df, dummy_encoded)
    F2 = AtomGraph(Float32.([0 1; 1 0]), ["F", "F"])
    @test_throws ArgumentError featurize!(F2, df)
end