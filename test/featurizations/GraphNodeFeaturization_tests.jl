using Test
using DataFrames
using CSV
using ChemistryFeaturization.FeatureDescriptor
using ChemistryFeaturization.Featurization

@testset "GraphNodeFeaturization" begin

    # custom lookup table...
    test_df = CSV.read(abspath(@__DIR__, "..", "test_data", "lookup_table.csv"), DataFrame)

    local fnames = ["X", "Block", "Atomic mass"]

    @testset "Encode" begin
        # make sure both constructors give the same results
        features = ElementFeatureDescriptor.(fnames)
        fzn1 = GraphNodeFeaturization(features)
        fzn2 = GraphNodeFeaturization(fnames)

        @test output_shape(fzn1) == 24 # 10 (default length for X) + 4 (default for block) + 10 (default for atomic mass)

        triangle_C_1 = AtomGraph(Float32.([0 1 1; 1 0 1; 1 1 0]), ["C", "C", "C"])
        triangle_C_2 = AtomGraph(Float32.([0 1 1; 1 0 1; 1 1 0]), ["C", "C", "C"])
        # @test fzn1 == fzn2

        featurized_1, featurized_2 = encode.([fzn1], [triangle_C_1, triangle_C_2]) # encode can be broadcasted
        @test featurized_1 == featurized_2

        featurized_1, featurized_2 = encode.([fzn1, fzn2], [triangle_C_1, triangle_C_2]) # encode can be broadcasted
        @test featurized_1 == featurized_2
    end

    @testset "Decode" begin
        # test that other options work properly

        # TODO - Figure out a way to robustly test this. Would creating a custom FD type whose encoded/decoded value you know beforehand and directly compare with be the best thing to do here?
        # fzn3 = GraphNodeFeaturization(fnames, nbins = 2)
        # F2 = AtomGraph(Float32.([0 1; 1 0]), ["F", "F"])

        # encode(F2, fzn3)
        # decoded_matrix = decode(fzn3, F2.encoded_features)
        # decoded_ag = decode(F2)
        # enc1 = F2.encoded_features
        # @test all(
        #     map(d -> d[1]["Block"] == d[2]["Block"] == "p", [decoded_matrix, decoded_ag]),
        # )
        # fzn4 = GraphNodeFeaturization(fnames, nbins = [2, 4, 2])
        # F2 = AtomGraph(Float32.([0 1; 1 0]), ["F", "F"])
        # encode(F2, fzn4)
        # @test all(F2.encoded_features .== enc1)
    end

    # encodable_elements
    @testset "Encodable Elements" begin
        fzn3 = GraphNodeFeaturization(ElementFeatureDescriptor.(["Boiling point", "6d"]))
        fzn4 = GraphNodeFeaturization(ElementFeatureDescriptor.(["7s", "Valence"]))
        @test encodable_elements(fzn3) == ["Ac", "Th", "U"]
        @test encodable_elements(fzn4) == ["Fr", "Ra", "Ac", "Th"]

        # Custom lookup_table
        feature_1 = ElementFeatureDescriptor("MeaningOfLife", test_df)   # zero-value case - `As` has a value = 0
        feature_2 = ElementFeatureDescriptor("neg_nums", test_df)
        feature_3 = ElementFeatureDescriptor("first_letter", test_df)
        feature_4 = ElementFeatureDescriptor("noarsenic", test_df) # missing value case - `As` has a missing value

        @test encodable_elements(GraphNodeFeaturization([feature_1, feature_2])) ==
              ["C", "As", "Tc"]
        @test encodable_elements(GraphNodeFeaturization([feature_1, feature_4])) ==
              ["C", "Tc"]
        @test encodable_elements(GraphNodeFeaturization([feature_2, feature_3])) ==
              ["C", "As", "Tc"]
        @test encodable_elements(GraphNodeFeaturization([feature_2, feature_4])) ==
              ["C", "Tc"]
    end

    # chunk_vec helper fcn
    @testset "chunk_vec" begin
        vec = [1, 1, 0, 1, 0, 1, 0]
        @test_throws AssertionError chunk_vec(vec, [3, 3])
        @test chunk_vec(vec, [4, 1, 2]) == [[1, 1, 0, 1], [0], [1, 0]]
    end

end
