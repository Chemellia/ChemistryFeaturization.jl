using Test
using DataFrames
using ..ChemistryFeaturization.Atoms: elements
using ..ChemistryFeaturization.FeatureDescriptor: OrbitalFeatureDescriptor
using SparseArrays

@testset "Encode-Decode" begin
    ofd = OrbitalFeatureDescriptor()

    @testset "Encode - with `default_ofd_encode`" begin
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

    @testset "Decode - with `default_ofd_decode`" begin
        alumina = AtomGraph(
            Float32.([0 1 0 0 0; 1 0 1 0 0; 0 1 0 1 0; 0 0 1 0 1; 0 0 0 1 0]),
            ["O", "Al", "O", "Al", "O"],
        )
        @test all(decode(ofd, ofd(alumina)) .== elements(alumina))
    end
end
