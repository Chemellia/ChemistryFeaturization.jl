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
    for i in 1:3
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
    df = DataFrame(:Symbol=>["C", "As"],:MeaningOfLife=>[42, 0])
    meaning = AtomFeature("MeaningOfLife", lookup_table=df)
    @test meaning(triangle_C)[10,:] == ones(3)
    @test encodable_elements(meaning) == ["C", "As"]

    # make a totally custom one, with a silly encode/decode
    element_encoder(element::String) = length(element)
    atoms_encoder = atoms -> reduce(hcat, map(e -> element_encoder(e), atoms.elements))
    decoder(encoded::Int) = string(['a' for i = 1:encoded]...)
    sillyfeature = AtomFeature("silly", atoms_encoder, decoder, false, false, 1, ["C", "Ne"])
    @test sillyfeature(triangle_C) == ones(1,3)
    @test encodable_elements(sillyfeature) == ["C", "Ne"]
end

#     # test that both constructors for categorical variables give the same result
#     cat1 = AtomFeat(:feat, true, 3, false, ['a','b','c'])
#     cat2 = AtomFeat(:feat, ['a','b','c'])
#     for name in propertynames(cat1)
#         @test getproperty(cat1,name)==getproperty(cat2,name)
#     end

#     # and now with categorical features that are numbers
#     cat3 = AtomFeat(:feat, true, 3, 1, 3)
#     cat3b = AtomFeat(:feat, true, 3, 1.0, 3)
#     cat4 = AtomFeat(:feat, true, 3, false, [1,2,3])
#     cat5 = AtomFeat(:feat, [1,2,3])

#     for name in [:name, :categorical, :num_bins, :logspaced]
#         results = getproperty.([cat3, cat3b, cat4, cat5], name)
#         @test all(y->y==results[1], results)
#     end

#     @test cat3.vals==cat4.vals==cat5.vals==[1,2,3]
#     @test cat3b.vals==[1.0,2.0,3.0]

#     # and for numerical features, first linearly spaced...
#     lin1 = AtomFeat(:feat, false, 3, false, [1.0,2,3,4])
#     lin2 = AtomFeat(:feat, false, 3, 1.0, 4)
#     @test lin1.vals==lin2.vals

#     # now log...
#     log1 = AtomFeat(:feat, false, 3, true, Float32.([0.1, 1, 10, 100]))
#     log2 = AtomFeat(:feat, false, 3, 0.1, 100, true)
#     @test log1.vals==log2.vals

#     # check that arbitrary spacing does work, though I don't see a need for it right now...
#     arb = AtomFeat(:feat, false, 3, false, [1,2,4,5])
#     @test arb.vals==Float32.([1,2,4,5])

#     # test build_atom_feats function
#     X, MP, block = build_featurization([:X, Symbol("Melting point"), :Block]; logspaced=[false, true])
#     @test X.name==:X
#     @test X.categorical==false
#     @test X.logspaced==false
#     @test X.vals[1]==Float32(0.7)
#     @test X.vals[end]==Float32(3.98)

#     @test MP.categorical==false
#     @test MP.logspaced==true

#     @test block.categorical==true
#     @test block.logspaced==false
#     @test block.vals==["s", "p", "d", "f"]
# end


# @testset "encode/decode" begin
#     # make_feature_vectors...pick some representative properties
#     feature_names = [:X, Symbol("Atomic mass"), :Block]
#     vecs, features = make_feature_vectors(feature_names)
#     # test that they are what we think they should be for a couple elements
#     H_feat = decode_feature_vector(vecs["H"], features)
#     @test H_feat[:X][1] <= 2.2 <= H_feat[:X][2]
#     @test H_feat[:Block] == "s"
#     @test H_feat[Symbol("Atomic mass")][1] <= Float32(1.00794) <= H_feat[Symbol("Atomic mass")][2]
#     Si_feat = decode_feature_vector(vecs["Si"], features)
#     @test Si_feat[:X][1] <= 1.9 <= Si_feat[:X][2]
#     @test Si_feat[:Block] == "p"
#     @test Si_feat[Symbol("Atomic mass")][1] <= 28.0855 <= Si_feat[Symbol("Atomic mass")][2]

#     # chunk_vec
#     vec = [1,0,0,1,0]
#     @test_throws AssertionError cf.chunk_vec(vec, [3,3])
#     @test cf.chunk_vec(vec, [3,2])==[[1,0,0],[1,0]]

#     # vec_valid
#     @test cf.vec_valid(vec, [3,2])
#     @test !cf.vec_valid(vec, [3,3])
#     @test !cf.vec_valid([1,0,1,1,0], [3,2])
# end
