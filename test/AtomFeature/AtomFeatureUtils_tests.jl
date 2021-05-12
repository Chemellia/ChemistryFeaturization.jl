using Test
const cf = ChemistryFeaturization

# @testset "binning" begin
# # test the which_bin function, can use a fake feature for numerical cases
# bins = [0,2,4]
# block_feat = AtomFeat(:Block, ["s", "p", "d", "f"])
# dummy_num_feat = AtomFeat(:Dummy, false, 2, false, [0,2,4])

# # test onehot_bins function
# @test cf.onehot_bins(block_feat, "s")==[1.0, 0., 0., 0.]
# @test cf.onehot_bins(dummy_num_feat, 1)==[1.0, 0.0]

# # and onecold_bins
# @test cf.onecold_bins(block_feat, [1,0,0,0])=="s"
# @test cf.onecold_bins(dummy_num_feat, [1, 0])==(0,2)

# # get_logspaced_vec
# @test cf.get_logspaced_vec(true, 3)==[true, true, true]
# @test cf.get_logspaced_vec(false, 3)==[false, false, false]
# @test cf.get_logspaced_vec([true,false,true], 3)==[true,false,true]
# @test cf.get_logspaced_vec([true,false,true], 2)==[true,false]
# end