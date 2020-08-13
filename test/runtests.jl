using CrystalGraphConvNets
using Test

@testset "SMILES_functions" begin
    include("SMILES_tests.jl")
    include("CIF_tests.jl")
end
