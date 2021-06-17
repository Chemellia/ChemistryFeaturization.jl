using Test
using DataFrames
using CSV
using ChemistryFeaturization.Utils.ElementFeatureUtils

@testset "ElementFeatureUtils" begin
    df = CSV.read(abspath(@__DIR__, "..", "test_data", "lookup_table.csv"), DataFrame)

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
    @test get_bins("Block") == ["s", "p", "d", "f"]
    @test get_bins("first_letter", df) == ["A", "C", "T"]
    @test get_bins("Group", categorical = false, nbins = 2) == 1.0:8.5:18.0
    @test get_bins("MeaningOfLife", df, categorical = true) == [-1, 0, 42]
    logbins = get_bins("Atomic mass")
    @test isapprox(logbins[2] / logbins[1], logbins[end] / logbins[end-1])
    neglogbins = get_bins("neg_nums", df, logspaced = true)
    @test neglogbins[1] == -1000
    # errors and things that should be ignored...
    @test get_bins("Block", nbins = 3) == ["s", "p", "d", "f"]
    @test get_bins("Valence", logspaced = true) == 0:1:13
    @test_throws AssertionError get_bins("Valence", categorical = false, logspaced = true)
    @test get_bins("noarsenic", df) == 1.0:0.1:2.0
    @test get_bins("noarsenic", df, categorical = true) == [1, 2]

    # build_onehot_vec
    @test build_onehot_vec("s", ["s", "p", "d", "f"], true) == [1, 0, 0, 0]
    @test_throws AssertionError build_onehot_vec("s", ["s", "p", "d", "f"], false)
    bins = 0:2:6
    @test build_onehot_vec(3, bins, false) == [0, 1, 0]
    @test build_onehot_vec(0, bins, false) == [1, 0, 0]
    @test build_onehot_vec(2, bins, false) == [0, 1, 0]
    @test build_onehot_vec(6, bins, false) == [0, 0, 1]
    @test_throws AssertionError build_onehot_vec(-1, bins, false)
    @test_throws ArgumentError build_onehot_vec(1, bins, true)

    # onehot_lookup_encoder
    @test onehot_lookup_encoder("C", "Block") == [0, 1, 0, 0]
    @test onehot_lookup_encoder("C", "Row")[2] == 1
    @test onehot_lookup_encoder("C", "MeaningOfLife", df, nbins = 4) == [0, 0, 0, 1]
    @test_throws TypeError onehot_lookup_encoder("He", "X")
    @test_throws TypeError onehot_lookup_encoder("As", "noarsenic", df)
    @test_throws AssertionError onehot_lookup_encoder("C", "heffalump")
    @test_throws AssertionError onehot_lookup_encoder("Si", "MeaningOfLife", df)

    # onecold_decoder
    @test onecold_decoder([0, 1, 0, 0], "Block") == "p"
    @test onecold_decoder([0 1; 1 0; 0 0; 0 0], "Block") == ["p", "s"]
    decoded = onecold_decoder(
        onehot_lookup_encoder("As", "MeaningOfLife", df),
        "MeaningOfLife",
        df,
    )
    @test decoded[1] <= 0 < decoded[2]
    @test_throws AssertionError onecold_decoder([1, 0, 0], "heffalump")

end
