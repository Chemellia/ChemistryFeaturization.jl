@testset "SimpleCodec" begin
    sc = SimpleCodec(x -> x .^ 2, x -> sqrt.(x))
    @test encode(3, sc) == 9
    @test encode([1, 2], sc) == [1, 4]
    @test decode(100, sc) == 10
    @test decode([1, 9], sc) == [1, 3]
end
