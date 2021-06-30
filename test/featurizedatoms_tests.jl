using Test: Serialization
using Test
using LightGraphs
using SimpleWeightedGraphs
using Serialization
using ..ChemistryFeaturization: FeaturizedAtoms


@testset "FeaturizedAtoms" begin

    local fnames = ["X", "Block", "Atomic mass"]

    # TODO - Actually test the `featurized_atoms.encoded_features` value. Probably define a custom type for this?
    @testset "Construct and encode" begin

        efds = ElementFeatureDescriptor.(fnames)
        fzn = GraphNodeFeaturization(efds)

        triangle_C = AtomGraph(Float32.([0 1 1; 1 0 1; 1 1 0]), ["C", "C", "C"])

        featurized_atoms = FeaturizedAtoms(triangle_C, fzn)

        # @test sef.encoded_features == ???
    end

    @testset "Decode" begin
        fzn = GraphNodeFeaturization(fnames, nbins = 2)
        F2 = AtomGraph(Float32.([0 1; 1 0]), ["F", "F"])
        featurized_atoms = FeaturizedAtoms(F2, fzn)

        decoded_matrix = decode(fzn, featurized_atoms.encoded_features)

        decoded_ag = decode(featurized_atoms)
        @test all(
            map(d -> d[1]["Block"] == d[2]["Block"] == "p", [decoded_matrix, decoded_ag]),
        )

        # check if `encoded_features` are generated correctly using encode itself
        fzn2 = GraphNodeFeaturization(fnames, nbins = [2, 4, 2])
        @test all(encode(fzn2, F2) .== featurized_atoms.encoded_features)

    end

end
