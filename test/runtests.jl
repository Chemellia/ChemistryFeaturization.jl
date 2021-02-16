using ChemistryFeaturization
using Test

@testset "weave_functions" begin
    include("weave_tests.jl")
end

@testset "Featurization" begin
    include("atomfeat_tests.jl")
end

@testset "Atomic Graphs" begin
    include("atomgraph_tests.jl")
end
