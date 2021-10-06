@testset "DirectCodec" begin
    dc = DirectCodec(1.0)
    @test encode(dc, 3) == 3
    @test encode(dc, [-2, 0]) == [-2, 0]
    @test decode(dc, 5) == 5
    @test decode(dc, [1, 2]) == [1, 2]
end
