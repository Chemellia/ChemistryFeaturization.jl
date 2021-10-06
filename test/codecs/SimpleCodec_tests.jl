@testset "SimpleCodec" begin
    sc = SimpleCodec(x -> x .^ 2, x -> sqrt.(x))
    @test encode(sc, 3) == 9
    @test encode(sc, [1, 2]) == [1, 4]
    @test decode(sc, 100) == 10
    @test decode(sc, [1, 9]) == [1, 3]
end
