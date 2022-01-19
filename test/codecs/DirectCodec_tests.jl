@testset "DirectCodec" begin
    dc = DirectCodec(1.0)
    @test encode(3, dc) == 3
    @test encode([-2, 0], dc) == [-2, 0]
    @test decode(5, dc) == 5
    @test decode([1, 2], dc) == [1, 2]
end
