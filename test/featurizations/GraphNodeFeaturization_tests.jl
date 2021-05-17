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

    # chunk_vec helper fcn

    # featurize!

    # decode

end