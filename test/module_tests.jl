using Test
using ChemistryFeaturization.AbstractType: AbstractFeaturization, AbstractFeatureDescriptor

@testset "modules and abstract methods" begin
    F2 = AtomGraph(Float32.([0 1; 1 0]), ["F", "F"])

    @testset "top-level module" begin
        # testing on `nothing` as example of ::Any
        @test_throws MethodError encodable_elements(nothing)
        @test_throws MethodError decode(nothing, nothing)
    end

    @testset "featurizations module" begin
        struct FakeFeaturization <: AbstractFeaturization end
        ff = FakeFeaturization()
        @test_throws MethodError encodable_elements(ff)
        @test_throws MethodError featurize!(F2, ff)
        @test_throws MethodError decode(ff, nothing)
    end

    # @testset "atoms module" begin
    # TBD cleanest way to test generic decode(::AbstractAtoms) - either another "fake" class, or maybe the `invoke` function
    # end

    @testset "features module" begin
        struct FakeFD <: AbstractFeatureDescriptor end
        fd = FakeFD()
        @test_throws MethodError encodable_elements(fd)
        @test_throws MethodError decode(fd, nothing)
    end
end
