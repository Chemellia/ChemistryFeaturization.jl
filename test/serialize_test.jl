using Test: Serialization
using Test
using LightGraphs
using SimpleWeightedGraphs
using Serialization
using ..ChemistryFeaturization: SerializableEncodedFeature


@testset "SerializableEncodedFeature" begin

    local fnames = ["X", "Block", "Atomic mass"]

    # TODO - Actually test the `sef.encoded_features` value. Probably define a custom type for this?
    @testset "Construct and featurize" begin

        efds = ElementFeatureDescriptor.(fnames)
        fzn = GraphNodeFeaturization(efds)

        triangle_C = AtomGraph(Float32.([0 1 1; 1 0 1; 1 1 0]), ["C", "C", "C"])

        sef = SerializableEncodedFeature(triangle_C, fzn)

        # @test sef.encoded_features == ???
    end

    @testset "Decode" begin
        fzn = GraphNodeFeaturization(fnames, nbins = 2)
        F2 = AtomGraph(Float32.([0 1; 1 0]), ["F", "F"])
        sef = SerializableEncodedFeature(F2, fzn)
    
        decoded_matrix = decode(fzn, sef.encoded_features)

        decoded_ag = decode(sef)
        @test all(
            map(d -> d[1]["Block"] == d[2]["Block"] == "p", [decoded_matrix, decoded_ag]),
        )

        # check if `encoded_features` are generated correctly using featurize itself
        fzn2 = GraphNodeFeaturization(fnames, nbins = [2, 4, 2])
        @test all(featurize(F2, fzn2) .== sef.encoded_features)

    end

end
