using Test
using DataFrames
using CSV

@testset "GraphNodeFeaturization" begin

    # make sure both constructors give the same results
    local fnames = ["X", "Block", "Atomic mass"]
    features = AtomFeature.(fnames)
    fzn1 = GraphNodeFeaturization(features)
    fzn2 = GraphNodeFeaturization(fnames)

    triangle_C_1 = AtomGraph(Float32.([0 1 1; 1 0 1; 1 1 0]), ["C", "C", "C"])
    triangle_C_2 = AtomGraph(Float32.([0 1 1; 1 0 1; 1 1 0]), ["C", "C", "C"])

    featurize!.([triangle_C_1, triangle_C_2], [fzn1, fzn2])
    @test triangle_C_1.atom_features == triangle_C_2.atom_features

    # test other options, make sure default values get appropriately populated...
    df = CSV.read(abspath(@__DIR__, "..", "test_data", "lookup_table.csv"), DataFrame)


    # encodable_elements
    @testset "Encodable Elements" begin
        fzn3 = GraphNodeFeaturization(AtomFeature.(["Boiling point", "6d"]))
        fzn4 = GraphNodeFeaturization(AtomFeature.(["7s", "Valence"]))
        @test encodable_elements(fzn3) == ["Ac", "Th", "U"]
        @test encodable_elements(fzn4) == ["Fr", "Ra", "Ac", "Th"]
    end


    # chunk_vec helper fcn
    @testset "chunk_vec" begin
        vec = [1, 1, 0, 1, 0, 1, 0]
        @test_throws AssertionError chunk_vec(vec, [3, 3])
        @test chunk_vec(vec, [4, 1, 2]) == [[1, 1, 0, 1], [0], [1, 0]]
    end

    # featurize!

    # decode

end
