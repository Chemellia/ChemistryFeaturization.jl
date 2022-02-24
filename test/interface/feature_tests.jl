struct DummyFeature <: AbstractFeatureDescriptor end

f = DummyFeature()

@testset "Features" begin
    @test_throws MethodError get_value(f, C3)
    @test_throws MethodError encodable_elements(f)
    @test_throws MethodError default_codec(f)
    @test_throws MethodError encode(C3, f)
end
