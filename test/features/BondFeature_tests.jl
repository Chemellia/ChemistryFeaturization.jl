using MolecularGraph

@testset "BondFeatureDescriptor" begin
    caffeine = smilestomol("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    ag = AtomGraph(caffeine)

    bfd = BondFeatureDescriptor("bondorder")
    vals = get_value(bfd, ag)
    @test vals[1,2] == vals[2,1] == 1
    @test vals[3,4] == vals[4,3] == 2
    @test all(ismissing.(vals[10,1:8]))
    # and encoding

    bfd = BondFeatureDescriptor("isringbond")
    vals = get_value(bfd, ag)
    @test vals[2,3] == vals[3,2] == true
    @test vals[14,9] == vals[9,14] == false
    # and encoding

    # build from scratch
end