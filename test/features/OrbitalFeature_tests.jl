using Test
using DataFrames
using ..ChemistryFeaturization: elements
using ..ChemistryFeaturization.FeatureDescriptor: OrbitalFeatureDescriptor
using SparseArrays

@testset "OrbitalFeatureDescriptor" begin
    ofd = OrbitalFeatureDescriptor()

    @testset "get_value" begin
        # individual elements
        @test get_value(ofd, "He") == "1s2"
        @test get_value(ofd, "C") == "2s2.2p2"
        @test get_value(ofd, "La") == "5d1.6s2"
        @test get_value(ofd, "Po") == "4f14.5d10.6s2.6p4"

        # and on an AtomGraph
        He_mol = AtomGraph(Float32.([1 0; 0 1]), ["He", "He"])
        @test get_value(ofd, He_mol) == ["1s2", "1s2"]
    end

    @testset "Encode" begin
        silica = AtomGraph(Float32.([0 1 0; 1 0 1; 0 1 0]), ["O", "Si", "O"])
        encoded_silica = encode(ofd, silica)
        I, J, V = findnz(encoded_silica)

        @test I == [2, 3, 4, 5, 2, 3]
        @test J == [1, 1, 2, 2, 3, 3]
        @test V == [2, 4, 2, 2, 2, 4]

        # non existent element - doesn't really make sense
        nonexistent = AtomGraph(Float32.([1 1; 1 1]), ["Non", "Existent"])
        @test_throws AssertionError ofd(nonexistent)
    end

    @testset "Decode" begin
        alumina = AtomGraph(
            Float32.([0 1 0 0 0; 1 0 1 0 0; 0 1 0 1 0; 0 0 1 0 1; 0 0 0 1 0]),
            ["O", "Al", "O", "Al", "O"],
        )
        @test all(decode(ofd, encode(ofd, alumina)) .== elements(alumina))
    end
end
