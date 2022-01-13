@testset "OneHotOneCold" begin
    ohoc_cat = OneHotOneCold(true, ['a', 'b', 'c'])
    ohoc_cont = OneHotOneCold(false, 0:3)
    @test_throws AssertionError OneHotOneCold(false, ['a', 'b', 'c'])

    @testset "Encode" begin
        @test encode(ohoc_cat, 'a') == [1, 0, 0]
        @test encode(ohoc_cat, ['a', 'c']) == [[1, 0, 0], [0, 0, 1]]
        @test_throws ArgumentError encode(ohoc_cat, 'd')

        @test encode(ohoc_cont, 0) == [1, 0, 0]
        @test encode(ohoc_cont, 1.5) == [0, 1, 0]
        @test encode(ohoc_cont, 3) == [0, 0, 1]
        @test_throws AssertionError encode(ohoc_cont, 4)
    end

    @testset "Decode" begin
        @test decode(ohoc_cat, [1, 0, 0]) == 'a'
        @test_throws AssertionError decode(ohoc_cat, [1, 0])

        @test decode(ohoc_cont, [0, 1, 0]) == (1, 2)
        @test decode(ohoc_cat, [1 1 1; 0 0 0; 0 0 0]) == ['a', 'a', 'a']
        @test decode(ohoc_cont, [1 1 1; 0 0 0; 0 0 0]) == [(0, 1), (0, 1), (0, 1)]
        arr = zeros(2, 2, 3)
        arr[:, :, 2] .= 1
        @test decode(ohoc_cat, arr) == ['b' 'b'; 'b' 'b']
        @test decode(ohoc_cont, arr) == [(1, 2) (1, 2); (1, 2) (1, 2)]
    end

end
