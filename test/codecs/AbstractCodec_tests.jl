struct DummyCodec <: AbstractCodec end

dc = DummyCodec()

@testset "AbstractCodec" begin
    @test_throws MethodError encode(0, dc)
    @test_throws MethodError decode(0, dc)
    @test_throws MethodError output_shape(dc)
    @test_throws MethodError output_shape(dc, 0)
end
