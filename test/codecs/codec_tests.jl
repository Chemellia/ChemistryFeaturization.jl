@testset "Codecs" begin
    codec_tests = [
        "AbstractCodec_tests",
        "OneHotOneCold_tests",
        "SimpleCodec_tests",
        "DirectCodec_tests",
    ]
    for t in codec_tests
        tp = abspath(testdir, "codecs", "$(t).jl")
        include(tp)
    end
end
