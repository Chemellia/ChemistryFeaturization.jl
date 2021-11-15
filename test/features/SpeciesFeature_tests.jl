using MolecularGraph
using ..ChemistryFeaturization.Codec

@testset "SpeciesFeatureDescriptor" begin
    # because if it doesn't work on caffeine then why bother...
    caffeine = smilestomol("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    ag = AtomGraph(caffeine)

    # because tea smells good...
    sfd = SpeciesFeatureDescriptor("isaromatic")
    @test all(get_value(sfd, ag) .== isaromatic(caffeine))
    @test length(encodable_elements(sfd)) == 118

    # test the other codec case...
    sfd = SpeciesFeatureDescriptor("lonepair")
    @test all(get_value(sfd, ag) .== lonepair(caffeine))
    encoded = encode(sfd, ag)
    @test all(encoded[1,:] .== encoded[5,:] .== 0)
    @test all(encoded[:, 8] .== encoded[:, 11] .== [0, 0, 0, 1, 0])

    # build one "from scratch"
    ag = AtomGraph(Float32.([0 1; 1 1]), ["H", "O"])
    num_nbs = g -> first.(length.(neighbors.(Ref(g), 1:nv(g))))
    codec = OneHotOneCold(true, [1, 2, 3, 4])
    categorical = true
    ee = ["O", "H"]
    sfd = SpeciesFeatureDescriptor{SimpleWeightedGraph,OneHotOneCold}(
        "num_neighbors",
        num_nbs,
        codec,
        categorical,
        ee,
    )
    @test get_value(sfd, ag) == [1, 2]
    @test encode(sfd, ag) == [1 0; 0 1; 0 0; 0 0]
    # add one to verify encodable_elements check happens
end
