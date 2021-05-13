using Test
using DataFrames
const cf = ChemistryFeaturization
const afu = ChemistryFeaturization.Utils.AtomFeatureUtils

@testset "AtomFeatureUtils" begin
    # one function at a time...

    # first make a dummy lookup table we'll use a few times...
    df = DataFrame(
        :Symbol => ["C", "As", "Tc"],
        :MeaningOfLife => [42, 0, -1],
        :first_letter => ['C', 'A', 'T'],
        :neg_nums => [-1, -10, -1000],
        :noarsenic => [1, missing, 2],
    )

    # fea_minmax
    @test_throws AssertionError fea_minmax("heffalump")
    @test fea_minmax("Group") == [1, 18]
    @test fea_minmax("MeaningOfLife", df) == [-1, 42]

    # default_log
    @test default_log("Block") == false # not numbers
    @test default_log("MeaningOfLife", df) == false # values span 0
    @test default_log("Valence") == false # extremal value 0
    @test default_log("Atomic mass") == true 
    @test default_log("Atomic mass", threshold = 3) == false
    # default_log
    @test default_log("Block") == false # not numbers
    @test default_log("MeaningOfLife", df) == false # values span 0
    @test default_log("Valence") == false # extremal value 0
    @test default_log("Atomic mass") == true
    @test default_log("Atomic mass", threshold = 3) == false

    # default_categorical
    @test default_categorical("Block") == true
    @test afu.get_param_vec(false, 2) == [false, false]
    @test afu.get_param_vec([false, false], 2) == [false, false]
    @test afu.get_param_vec([false, false], 3) == [false, false, false]
    @test afu.get_param_vec([false, false, false], 2) == [false, false]
    @test afu.get_param_vec([3, 2], 3, pad_val = 0) == [3, 2, 0]
    @test afu.get_param_vec(['a'], 3, pad_val = 'b') == ['a', 'b', 'b']
    @test afu.get_param_vec('a', 3, pad_val = 'b') == ['a','a','a']

    @test afu.get_param_vec([false, false], 3) == [false, false, false]
    @test afu.get_param_vec([false, false, false], 2) == [false, false]
    @test afu.get_param_vec([3, 2], 3, pad_val = 0) == [3, 2, 0]
    @test afu.get_param_vec(['a'], 3, pad_val = 'b') == ['a', 'b', 'b']
    @test afu.get_param_vec('a', 3, pad_val = 'b') == ['a', 'a', 'a']

    # get_bins
    @test get_bins("Block") == ["s", "p", "d", "f"]
    @test get_bins("first_letter", df) == ['A', 'C', 'T']
    @test get_bins("Group", categorical = false, nbins = 2) == 1.0:8.5:18.0
    @test get_bins("MeaningOfLife", df, categorical = true) == [-1, 0, 42]
    logbins = get_bins("Atomic mass")
    @test logbins[2] / logbins[1] == logbins[end] / logbins[end-1]
    neglogbins = get_bins("neg_nums", df, logspaced = true)
    @test neglogbins[1] == -1000
    # errors and things that should be ignored...
    @test get_bins("Block", nbins = 3) == ["s", "p", "d", "f"]
    @test get_bins("Valence", logspaced = true) == 0:1:13
    @test_throws AssertionError get_bins("Valence", categorical = false, logspaced = true)
    @test get_bins("noarsenic", df) == 1.0:0.1:2.0
    @test get_bins("noarsenic", df, categorical = true) == [1, 2]

    # build_onehot_vecs

    # onehot_lookup_encoder

    # onecold_decoder

    # encodable_elements

end