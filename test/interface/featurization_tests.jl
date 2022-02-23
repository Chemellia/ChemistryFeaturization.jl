@testset "Featurization" begin
    @test !("He" in encodable_elements(fzn1))
    @test_throws MethodError features(fzn2)

    encoded = encode(C3, fzn1)
    @test size(encoded[1]) == (4, 3)
    @test size(encoded[2]) == (10, 3)
    @test all(encoded[1][3, :] .== 1)

    decoded = decode(encoded, fzn1)
    @test all(decoded[1] .== "p")
    @test decoded[2][1] == (2.34, 2.668)
end
