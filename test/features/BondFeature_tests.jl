using MolecularGraph
using SimpleWeightedGraphs
using ..ChemistryFeaturization.Codec

@testset "BondFeatureDescriptor" begin
    caffeine = smilestomol("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
    ag = AtomGraph(caffeine)

    bfd = BondFeatureDescriptor("bondorder")
    vals = get_value(bfd, ag)
    @test vals[1, 2] == vals[2, 1] == 1
    @test vals[3, 4] == vals[4, 3] == 2
    @test all(ismissing.(vals[10, 1:8]))
    encoded = encode(bfd, ag)
    @test encoded[4, 5, 1] == encoded[5, 4, 1] == true
    @test all(encoded[4, 5, 2:3] .== encoded[5, 4, 2:3] .== false)
    @test all(skipmissing(decode(bfd, encoded) .== vals))

    bfd = BondFeatureDescriptor("isringbond")
    vals = get_value(bfd, ag)
    @test vals[2, 3] == vals[3, 2] == true
    @test vals[14, 9] == vals[9, 14] == false
    @test encode(bfd, ag) .== vals
    @test all(skipmissing(decode(bfd, encoded) .== vals))

    # build from scratch
    ag = AtomGraph(Float32.([0 1; 1 0.5]), ["H", "O"])
    function edge_wt(g::SimpleWeightedGraph)
        wts = collect(weights(g))
        output = Matrix{Union{Float64,Missing}}(missing, size(wts)...)
        for ind in CartesianIndices(wts)
            if !(wts[ind] == 0)
                output[ind] = wts[ind]
            end
        end
        return output
    end
    categorical = false
    codec = OneHotOneCold(categorical, [0, 1 / 3, 2 / 3, 1])
    ee = ["O", "H"]
    bfd = SpeciesFeatureDescriptor{SimpleWeightedGraph,OneHotOneCold}(
        "edge_weight",
        edge_wt,
        codec,
        categorical,
        ee,
    )
    vals = get_value(bfd, ag)
    @test vals[1, 2] == vals[2, 1] == 1.0
    @test ismissing(vals[1, 1])
    @test vals[2, 2] == 0.5
    # NB no test for decoding here right now beacuse there's no dispatch onto get_value of BFD for AtomGraph{SimpleWeightedGraph}
end
