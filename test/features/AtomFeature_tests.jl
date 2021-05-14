using Test
using DataFrames
const cf = ChemistryFeaturization

@testset "AtomFeature" begin
    # test for errors...
    @test_throws AssertionError AtomFeature("heffalump")

    # construct a few from built-in features...
    local fnames = ["X", "Block", "Atomic mass"]
    local lengths = [default_nbins, 4, default_nbins]
    local cats = [false, true, false]
    local conts = [false, false, false]
    for i = 1:3
        f = AtomFeature(fnames[i])
        @test f.length == lengths[i]
        @test f.categorical == cats[i]
        @test f.contextual == conts[i]
    end

    # encode/decode some stuff...
    He_mol = AtomGraph(Float32.([1 0; 0 1]), ["He", "He"])
    X, block, amass = AtomFeature.(fnames)
    @test_throws AssertionError X(He_mol)
    @test block(He_mol) == Float64.([0 0; 1 1; 0 0; 0 0])
    @test decode(block, block(He_mol)) == ["p", "p"]
    @test amass(He_mol)[3, :] == ones(2)
    true_He_amass = atom_data_df[2, Symbol("Atomic mass")]
    He_amass_min, He_amass_max = decode(amass, amass(He_mol)[:, 1])
    @test He_amass_min < true_He_amass < He_amass_max

    # now let's try some options (and make sure they're ignored when appropriate...)...
    X = AtomFeature("X", nbins = 8, logspaced = true)
    block = AtomFeature("Block", nbins = 5)
    triangle_C = AtomGraph(Float32.([0 1 1; 1 0 1; 1 1 0]), ["C", "C", "C"])
    @test X(triangle_C)[6, :] == ones(3)
    # TODO: decoded values
    @test block(triangle_C)[2, :] == ones(3)

    # and make a custom lookup table...
    df = DataFrame(:Symbol => ["C", "As"], :MeaningOfLife => [42, 0])
    meaning = AtomFeature("MeaningOfLife", df)
    @test meaning(triangle_C)[10, :] == ones(3)
    @test encodable_elements(meaning) == ["C", "As"]
    @test encodable_elements("MeaningOfLife", df) == ["C", "As"]

    # make a totally custom one, with a silly encode/decode
    element_encoder(element::String) = length(element)
    atoms_encoder = atoms -> reduce(hcat, map(e -> element_encoder(e), atoms.elements))
    decoder(encoded::Int) = string(['a' for i = 1:encoded]...)
    sillyfeature =
        AtomFeature("silly", atoms_encoder, decoder, false, false, 1, ["C", "Ne"])
    @test sillyfeature(triangle_C) == ones(1, 3)
    @test encodable_elements(sillyfeature) == ["C", "Ne"]
end
