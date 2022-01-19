using DataFrames
using CSV
using ChemistryFeaturization.ElementFeature

@testset "ElementFeatureUtils" begin
    df = CSV.read(abspath(@__DIR__, "test_data", "lookup_table.csv"), DataFrame)

    # fea_minmax
    @test_throws AssertionError ElementFeature.fea_minmax("heffalump")
    @test ElementFeature.fea_minmax("Group") == [1, 18]
    @test ElementFeature.fea_minmax("MeaningOfLife", df) == [-1, 42]

    # default_log
    @test default_log("Block") == false # not numbers
    @test default_log("MeaningOfLife", df) == false # values span 0
    @test default_log("Valence") == false # extremal value 0
    @test default_log("Atomic mass") == true
    @test default_log("Atomic mass", threshold_oom = 3) == false

    # default_categorical
    @test default_categorical("Block") == true
    @test get_param_vec(false, 2) == [false, false]
    @test get_param_vec([false, false], 2) == [false, false]
    @test get_param_vec([false, false], 3) == [false, false, false]
    @test get_param_vec([false, false, false], 2) == [false, false]
    @test get_param_vec([3, 2], 3, pad_val = 0) == [3, 2, 0]
    @test get_param_vec(['a'], 3, pad_val = 'b') == ['a', 'b', 'b']
    @test get_param_vec('a', 3, pad_val = 'b') == ['a', 'a', 'a']

    @test get_param_vec([false, false], 3) == [false, false, false]
    @test get_param_vec([false, false, false], 2) == [false, false]
    @test get_param_vec([3, 2], 3, pad_val = 0) == [3, 2, 0]
    @test get_param_vec(['a'], 3, pad_val = 'b') == ['a', 'b', 'b']
    @test get_param_vec('a', 3, pad_val = 'b') == ['a', 'a', 'a']

    # get_bins
    @test get_bins("Block") == ["d", "f", "p", "s"]
    @test get_bins("first_letter", df) == ["A", "C", "T"]
    @test get_bins("Group", categorical = false, nbins = 2) == 1.0:8.5:18.0
    @test get_bins("MeaningOfLife", df, categorical = true) == [-1, 0, 42]
    logbins = get_bins("Atomic mass")
    @test isapprox(logbins[2] / logbins[1], logbins[end] / logbins[end-1])
    neglogbins = get_bins("neg_nums", df, logspaced = true)
    @test neglogbins[1] == -1000
    # errors and things that should be ignored...
    @test get_bins("Block", nbins = 3) == ["d", "f", "p", "s"]
    @test get_bins("Valence", logspaced = true) == 0:1:13
    @test_throws AssertionError get_bins("Valence", categorical = false, logspaced = true)
    @test get_bins("noarsenic", df) == [1, 2]
    @test get_bins("noarsenic", df, categorical = true) == [1, 2]

end
