using Test
using DataFrames
using CSV

@testset "GraphNodeFeaturization" begin

    # make sure both constructors give the same results
    local fnames = ["X", "Block", "Atomic mass"]
    features = ElementFeatureDescriptor.(fnames)
    fzn1 = GraphNodeFeaturization(features)
    fzn2 = GraphNodeFeaturization(fnames)

    triangle_C_1 = AtomGraph(Float32.([0 1 1; 1 0 1; 1 1 0]), ["C", "C", "C"])
    triangle_C_2 = AtomGraph(Float32.([0 1 1; 1 0 1; 1 1 0]), ["C", "C", "C"])

    # test other options, make sure default values get appropriately populated...
    test_df = CSV.read(abspath(@__DIR__, "..", "test_data", "lookup_table.csv"), DataFrame)

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

    # featurize!
    featurize!.([triangle_C_1, triangle_C_2], [fzn1, fzn2])
    @test triangle_C_1.encoded_atom_features == triangle_C_2.encoded_atom_features

    # decode

end
