using DataFrames
using CSV
using ..ChemistryFeaturization.ElementFeature
using ..ChemistryFeaturization.Data

@testset "ElementFeatureDescriptor" begin
    # construct a few from built-in features...
    X = ElementFeatureDescriptor("X")
    custom_xc = OneHotOneCold(false, get_bins("X", nbins=8, logspaced=true))
    block = ElementFeatureDescriptor("Block")
    amass = ElementFeatureDescriptor("Atomic mass")
    @testset "default codec" begin
        xc = default_codec(X)
        @test !xc.categorical
        @test xc.bins[1] == 0.7
        @test xc.bins[end] == 3.98

        bc = default_codec(block)
        @test bc.categorical
        @test bc.bins == ["d", "f", "p", "s"]
    end

    @testset "Encode" begin
        @test get_value(block, C3) == block(C3)
        @test_throws AssertionError get_value(X, He_mol)

        @test_throws AssertionError encode(X, He_mol)
        @test encode(block, He_mol) == Float64.([0 0; 0 0; 1 1; 0 0])
        @test encode(amass, He_mol)[3, :] == ones(2)

        # explicit codec
        @test encode(X, custom_xc, C3)[6, :] == ones(3)
    end

    @testset "Decode" begin
        @test decode(block, encode(block, He_mol)) == ["p", "p"]

        true_He_amass = element_data_df[2, Symbol("Atomic mass")]
        He_amass_min, He_amass_max = decode(amass, encode(amass, He_mol)[:, 1])
        @test He_amass_min < true_He_amass < He_amass_max

        true_X = element_data_df.X[6]
        decoded_X_range = decode(X, encode(X, C3)[:, 1])
        @test decoded_X_range[1] < true_X < decoded_X_range[2]
        @test decode(block, encode(block, C3)) == ["p", "p", "p"]
    end

    # and make a custom lookup table...
    @testset "Custom Lookup Table" begin
        df = CSV.read(abspath(@__DIR__, "test_data", "lookup_table.csv"), DataFrame)
        meaning = ElementFeatureDescriptor("MeaningOfLife", df)
        @test encode(meaning, C3)[3, :] == ones(3)

        @test encodable_elements(meaning) == ["C", "As", "Tc"]
        @test_throws AssertionError get_value(meaning, He_mol)
    end

    @testset "errors" begin
        @test_throws AssertionError ElementFeatureDescriptor("heffalump")
    end
end
