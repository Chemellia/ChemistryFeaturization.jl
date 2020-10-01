using ChemistryFeaturization
using Test

@testset "SMILES_functions" begin
    include("SMILES_tests.jl")
end

@testset "Featurization" begin
    include("atomfeat_tests.jl")
end

@testset "Atomic Graphs" begin
    include("atomgraph_tests.jl")
end
