using Test

@testset "SMILES_functions" begin
    include("SMILES_tests.jl")
end

@testset "CIF_functions" begin
    include("CIF_tests.jl")
end

@testset "types" begin
    include("type_tests.jl")
end