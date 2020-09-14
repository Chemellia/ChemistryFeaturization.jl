using Test
const cf = ChemistryFeaturization

@testset "binning" begin
    # test the which_bin function, can use a fake feature for numerical cases
    bins = [0,2,4]
    block_feat = AtomFeat(:Block, ["s", "p", "d", "f"])
    dummy_num_feat = AtomFeat(:Dummy, false, 2, false, [0,2,4])

    # test onehot_bins function
    @test cf.onehot_bins(block_feat, "s")==[1.0, 0., 0., 0.]
    @test cf.onehot_bins(dummy_num_feat, 1)==[1.0, 0.0]

    # and onecold_bins
    @test cf.onecold_bins(block_feat, [1,0,0,0])=="s"
    @test cf.onecold_bins(dummy_num_feat, [1, 0])==(0,2)

    # get_logspaced_vec
    @test cf.get_logspaced_vec(true, 3)==[true, true, true]
    @test cf.get_logspaced_vec(false, 3)==[false, false, false]
    @test cf.get_logspaced_vec([true,false,true], 3)==[true,false,true]
    @test cf.get_logspaced_vec([true,false,true], 2)==[true,false]
end

@testset "AtomFeat" begin
    # checks that bins get created correctly
    # first for categorical features, of at least a couple different types

    # and for numerical features, both linear and log spaced

    # test all the constructors and functions obviously...
end

@testset "encode/decode" begin
    # make_feature_vectors...pick some representative properties
    feature_names = [:X, Symbol("Atomic mass"), :Block]
    vecs, features = make_feature_vectors(feature_names)
    # test that they are what we think they should be for a couple elements
    H_feat = decode_feature_vector(vecs["H"], features)
    @test H_feat[:X][1] <= 2.2 <= H_feat[:X][2]
    @test H_feat[:Block] == "s"
    @test H_feat[Symbol("Atomic mass")][1] <= 1.00794 <= H_feat[Symbol("Atomic mass")][2]
    Si_feat = decode_feature_vector(vecs["Si"], features)
    @test Si_feat[:X][1] <= 1.9 <= Si_feat[:X][2]
    @test Si_feat[:Block] == "p"
    @test Si_feat[Symbol("Atomic mass")][1] <= 28.0855 <= Si_feat[Symbol("Atomic mass")][2]

    # chunk_vec
    vec = [1,0,0,1,0]
    @test_throws AssertionError cf.chunk_vec(vec, [3,3])
    @test cf.chunk_vec(vec, [3,2])==[[1,0,0],[1,0]]

    # vec_valid
    @test cf.vec_valid(vec, [3,2])
    @test !cf.vec_valid(vec, [3,3])
    @test !cf.vec_valid([1,0,1,1,0], [3,2])
end