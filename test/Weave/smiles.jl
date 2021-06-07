using DataStructures
using Test
using ChemistryFeaturization
include("../../src/weave_fcns.jl")
using .weave_fcns:
    smiles_atom_features,
    smiles_bond_features,
    chem,
    atom_feat_list,
    bond_feat_bins,
    atom_feat_bins


mol = chem.MolFromSmiles("C1=CC=CC=C1")  # benzene
carbon = get(mol.GetAtoms(), 0)
bond = get(mol.GetBonds(), 0)


@testset "SMILES atom featurization" begin
    @test [k for k in keys(smiles_atom_features(carbon))] == [
        "symbol",
        "degree",
        "implicit_valence",
        "formal_charge",
        "radical_electrons",
        "hybridization",
        "aromaticity",
        "total_H_num",
    ]
    try
        smiles_atom_features(carbon, feature_list = ["invalid_feat"])
    catch e
        @test occursin("is not in supported atomic features list", sprint(showerror, e))
    end
    correct_values = [
        Float32.("C" .== atom_feat_bins["symbol"]),
        Float32.(2 .== atom_feat_bins["degree"]),
        Float32.(1 .== atom_feat_bins["implicit_valence"]),
        Float32(carbon.GetFormalCharge()),
        Float32(carbon.GetNumRadicalElectrons()),
        Float32.(chem.rdchem.HybridizationType.SP2 .== atom_feat_bins["hybridization"]),
        Float32(carbon.GetIsAromatic()),
        Float32.(1 .== atom_feat_bins["total_H_num"]),
    ]
    @test [v for v in values(smiles_atom_features(carbon))] == correct_values
end


@testset "SMILES bond featurization" begin
    @test [k for k in keys(smiles_bond_features(bond))] == ["bond_type", "isConjugated", "isInring"]
    try
        smiles_bond_features(carbon, feature_list = ["invalid_feat"])
    catch e
        @test occursin("is not in supported bond features list", sprint(showerror, e))
    end
    correct_values = [
        Float32.(chem.rdchem.BondType.AROMATIC .== bond_feat_bins["bond_type"]),
        Float32(bond.GetIsConjugated()),
        Float32(bond.IsInRing()),
    ]
    @test [v for v in values(smiles_bond_features(bond))] == correct_values
end
