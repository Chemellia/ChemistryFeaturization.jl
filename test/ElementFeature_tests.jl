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

        @test_throws AssertionError encode(He_mol, X)
        @test encode(He_mol, block) == Float64.([0 0; 0 0; 1 1; 0 0])
        @test encode(He_mol, amass)[3, :] == ones(2)

        # explicit codec
        @test encode(C3, X, custom_xc)[6, :] == ones(3)
    end

    @testset "Decode" begin
        @test decode(encode(He_mol, block), block) == ["p", "p"]

        true_He_amass = element_data_df[2, Symbol("Atomic mass")]
        He_amass_min, He_amass_max = decode(encode(He_mol, amass)[:, 1], amass)
        @test He_amass_min < true_He_amass < He_amass_max

        true_X = element_data_df.X[6]
        decoded_X_range = decode(encode(C3, X)[:, 1], X)
        @test decoded_X_range[1] < true_X < decoded_X_range[2]
        @test decode(encode(C3, block), block) == ["p", "p", "p"]
    end

    # and make a custom lookup table...
    @testset "Custom Lookup Table" begin
        df = CSV.read(abspath(@__DIR__, "test_data", "lookup_table.csv"), DataFrame)
        meaning = ElementFeatureDescriptor("MeaningOfLife", df)
        @test encode(C3, meaning)[3, :] == ones(3)

        @test encodable_elements(meaning) == ["C", "As", "Tc"]
        @test_throws AssertionError get_value(meaning, He_mol)
    end

    @testset "errors" begin
        @test_throws AssertionError ElementFeatureDescriptor("heffalump")
    end
end
