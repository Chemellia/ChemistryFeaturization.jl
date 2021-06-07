using DataStructures
using Test
using ChemistryFeaturization
include("../../src/featurize.jl")
using .featurize: weave_featurize

mol = chem.MolFromSmiles("C1=CC=CC=C1")  # benzene
carbon = get(mol.GetAtoms(), 0)
bond = get(mol.GetBonds(), 0)


@testset "weave featurization" begin
    weave_mol = weave_featurize(["C1=CC=CC=C1"])[1]
    @test sum(weave_mol.atom_features[end, :]) ==
          sum([sum(v) for v in values(smiles_atom_features(carbon))])
    correct_pair_feat_sum = sum([sum(v) for v in values(smiles_bond_features(bond))]) + 2
    count = 0
    for i = 1:size(weave_mol.pair_features, 1)
        if sum(weave_mol.pair_features[i, :]) == correct_pair_feat_sum
            count += 1
        end
    end
    @test count == 12
    @test weave_mol.num_atoms == 6
    @test weave_mol.atom_feat_num == 75

    @test_logs (:warn, "No bonds found in C! Pair feature matrix filled with zeros!") weave_featurize([
        "C",
    ])[1]

end
