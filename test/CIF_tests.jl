using Test

include("../src/pmg_features.jl")
include("../src/pmg_graphs.jl")

@testset "binning" begin
    # test the get_bins function
    # first for categorical features
    for feature in keys(categorical_feature_vals)
        @test get_bins(feature)==categorical_feature_vals[feature]
    end
    # TODO: add tests for numerical features

    # test the which_bin function, can use a fake feature for numerical cases
    bins = [0,2,4]
    @test which_bin("", 0, bins)==1
    @test which_bin("", 1, bins)==1
    @test which_bin("", 3, bins)==2
    @test which_bin("", 4, bins)==2
    # TODO: add test for categorical features

    # test onehot_bins function
    @test onehot_bins(:Block, "s")==[1.0, 0., 0., 0.]
    @test onehot_bins("", 1, bins)==[1.0, 0.0]

    # and onecold_bins
    @test onecold_bins(:Block, [1,0,0,0], ["s","p","d","f"])=="s"
    @test onecold_bins("", [1, 0], bins)==(0,2)

    # get_logspaced_vec
    @test get_logspaced_vec(true, 3)==[true, true, true]
    @test get_logspaced_vec(false, 3)==[false, false, false]
    @test get_logspaced_vec([true,false,true], 3)==[true,false,true]
    @test get_logspaced_vec([true,false,true], 2)==[true,false]
end

@testset "encode/decode" begin
    # make_feature_vectors...pick some representative properties
    features = [:X, Symbol("Atomic mass"), :Block]
    vecs = make_feature_vectors(features)
    # test that they are what we think they should be for a couple elements
    H_feat = decode_feature_vector(vecs["H"], features, [default_nbins, default_nbins, 4])
    @test H_feat[:X][1] <= 2.2 <= H_feat[:X][2]
    @test H_feat[:Block] == "s"
    @test H_feat[Symbol("Atomic mass")][1] <= 1.00794 <= H_feat[Symbol("Atomic mass")][2]
    Si_feat = decode_feature_vector(vecs["Si"], features, [default_nbins, default_nbins, 4])
    @test Si_feat[:X][1] <= 1.9 <= Si_feat[:X][2]
    @test Si_feat[:Block] == "p"
    @test Si_feat[Symbol("Atomic mass")][1] <= 28.0855 <= Si_feat[Symbol("Atomic mass")][2]

    # chunk_vec
    vec = [1,0,0,1,0]
    @test_throws AssertionError chunk_vec(vec, [3,3])
    @test chunk_vec(vec, [3,2])==[[1,0,0],[1,0]]

    # vec_valid
    @test vec_valid(vec, [3,2])
    @test !vec_valid(vec, [3,3])
    @test !vec_valid([1,0,1,1,0], [3,2])
end

@testset "graph-building" begin
wm, atoms = build_graph(joinpath(@__DIR__, "./test_data/mp-195.cif"))
    wm_true = [0.0 1.0 1.0 1.0; 1.0 0.0 1.0 1.0; 1.0 1.0 0.0 1.0; 1.0 1.0 1.0 0.0]
    @test wm == wm_true
    @test atoms == ["Ho", "Pt", "Pt", "Pt"]
    wm, atoms = build_graph(joinpath(@__DIR__, "./test_data/mp-195.cif"; use_voronoi=false)
    @test wm == wm_true
    @test atoms == ["Ho", "Pt", "Pt", "Pt"]
end
