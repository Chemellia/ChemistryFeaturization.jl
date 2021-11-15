using Test
using ChemistryFeaturization.AbstractType:
    AbstractAtoms, AbstractFeaturization, AbstractFeatureDescriptor, AbstractCodec
using ChemistryFeaturization.Atoms: elements
using ChemistryFeaturization: elements

@testset "Modules and Abstract methods" begin
    struct FakeAtoms <: AbstractAtoms{Nothing} end
    struct FakeFD <: AbstractFeatureDescriptor end
    struct FakeFeaturization <: AbstractFeaturization end
    struct FakeCodec <: AbstractCodec end

    @testset "top-level module" begin
        # testing on `nothing` as example of ::Any
        @test_throws MethodError encodable_elements(nothing)
        @test_throws MethodError encode(nothing, nothing)
        @test_throws MethodError decode(nothing, nothing)
        @test_throws MethodError elements(nothing)
    end

    @testset "atoms module" begin
        fa = FakeAtoms()
        @test_throws MethodError elements(fa)
    end

    @testset "features module" begin
        fa = FakeAtoms()
        ffd = FakeFD()
        fc = FakeCodec()
        @test_throws ErrorException encodable_elements(ffd)
        @test_throws MethodError get_value(ffd, fa)
        @test_throws ErrorException encode(ffd, fa)
        @test_throws ErrorException decode(ffd, nothing)
    end

    @testset "featurizations module" begin
        ff = FakeFeaturization()
        F2_atom = AtomGraph(Float32.([0 1; 1 0]), ["F", "F"])

        @test_throws MethodError encodable_elements(ff)
        @test_throws MethodError encode(ff, F2_atom)
        @test_throws MethodError decode(ff, nothing)
    end
end
