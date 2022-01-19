@testset "OneHotOneCold" begin
    ohoc_cat = OneHotOneCold(true, ['a', 'b', 'c'])
    ohoc_cont = OneHotOneCold(false, 0:3)
    @test_throws AssertionError OneHotOneCold(false, ['a', 'b', 'c'])

    @testset "Encode" begin
        @test encode('a', ohoc_cat) == [1, 0, 0]
        @test encode(['a', 'c'], ohoc_cat) == [[1, 0, 0], [0, 0, 1]]
        @test_throws ArgumentError encode('d', ohoc_cat)

        @test encode(0, ohoc_cont) == [1, 0, 0]
        @test encode(1.5, ohoc_cont) == [0, 1, 0]
        @test encode(3, ohoc_cont) == [0, 0, 1]
        @test_throws AssertionError encode(4, ohoc_cont)
    end

    @testset "Decode" begin
        @test decode([1, 0, 0], ohoc_cat) == 'a'
        @test_throws AssertionError decode([1, 0], ohoc_cat)

        @test decode([0, 1, 0], ohoc_cont) == (1, 2)
        @test decode([1 1 1; 0 0 0; 0 0 0], ohoc_cat) == ['a', 'a', 'a']
        @test decode([1 1 1; 0 0 0; 0 0 0], ohoc_cont) == [(0, 1), (0, 1), (0, 1)]
        arr = zeros(2, 2, 3)
        arr[:, :, 2] .= 1
        @test decode(arr, ohoc_cat) == ['b' 'b'; 'b' 'b']
        @test decode(arr, ohoc_cont) == [(1, 2) (1, 2); (1, 2) (1, 2)]
    end

end
