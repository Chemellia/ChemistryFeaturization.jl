using Test
using DataFrames
using CSV
using ..ChemistryFeaturization.FeatureDescriptor
const cf = ChemistryFeaturization

@testset "ElementFeatureDescriptor" begin
    # test for errors
    @test_throws AssertionError ElementFeatureDescriptor("heffalump")

    # construct a few from built-in features...
    local fnames = ["X", "Block", "Atomic mass"]
    local lengths = [default_nbins, 4, default_nbins]
    local cats = [false, true, false]
    local conts = [false, false, false]

    for i = 1:3
        f = ElementFeatureDescriptor(fnames[i])
        @test output_shape(f) == lengths[i]
        @test f.categorical == cats[i]
    end

    # # encode/decode some stuff...
    He_mol = AtomGraph(Float32.([1 0; 0 1]), ["He", "He"])
    X, block, amass = ElementFeatureDescriptor.(fnames)
    @test_throws AssertionError X(He_mol)

    @test block(He_mol) == Float64.([0 0; 1 1; 0 0; 0 0])

    @test decode(block, block(He_mol)) == ["p", "p"]

    @test amass(He_mol)[3, :] == ones(2)

    true_He_amass = atom_data_df[2, Symbol("Atomic mass")]
    He_amass_min, He_amass_max = decode(amass, amass(He_mol)[:, 1])
    @test He_amass_min < true_He_amass < He_amass_max

    # now let's try some options (and make sure they're ignored when appropriate...)...
    X = ElementFeatureDescriptor("X", nbins = 8, logspaced = true)
    triangle_C = AtomGraph(Float32.([0 1 1; 1 0 1; 1 1 0]), ["C", "C", "C"])
    @test X(triangle_C)[6, :] == ones(3)

    block = ElementFeatureDescriptor("Block", nbins = 5)
    @test block(triangle_C)[2, :] == ones(3)

    true_X = atom_data_df.X[6]
    decoded_X_range = decode(X, X(triangle_C)[:, 1])
    @test decoded_X_range[1] < true_X < decoded_X_range[2]
    @test decode(block, block(triangle_C)) == ["p", "p", "p"]

    # and make a custom lookup table...
    df = CSV.read(abspath(@__DIR__, "..", "test_data", "lookup_table.csv"), DataFrame)
    meaning = ElementFeatureDescriptor("MeaningOfLife", df)
    @test meaning(triangle_C)[10, :] == ones(3)
    @test encodable_elements(meaning) == ["C", "As", "Tc"]
    @test encodable_elements("MeaningOfLife", df) == ["C", "As", "Tc"]
end
