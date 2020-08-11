# this file will have the "master" function that's actually exposed to the user
include("cif_fcns.jl")
include("smiles_fcns.jl")

function featurize_inputs(input_data_folder, SMILES_atom_features, SMILES_bond_features, other_atom_features, graphs_folder, features_folder)
    # determine type of data in input folder (must be all SMILES or all CIFs)
    # check whether user has asked for features that make sense (may require some ancillary fcns)
    # loop over input files and farm out actual parsing to specialized functions for CIF vs. SMILES
    # eventually, save out adjacency matrices to graphs_folder and feature vectors/metadata to features_folder
end


#downstream methods for concatnating atom features generate by using SMILES
concat_smiles_atom_feat(atom_feat::OrderedDict) = vcat(values(atom_feat)...)
concat_smiles_bond_feat(bond_feat::OrderedDict) = vcat(values(bond_feat)...)
## Weave Featurization

# data structure for storing features returned from wevae_featurizer
mutable struct WeaveMol
    atom_features
    pair_features
    num_atoms
    atom_feat_num
end


"""
weave featurizier

input:
    inputs(Array{String}): List of SIMLES strings
    atom_feature_list(Array{String}): list of atom features of interest
    bond_feature_list(Array{String}): list of bond features of interest
"""
function weave_featurize(inputs; atom_feature_list::Array{String}=atom_feat_list, bond_feature_list::Array{String}=bond_feat_list)
    outputs = Array{WeaveMol}(undef,length(inputs))
    for (i,smiles) in enumerate(inputs)
        mol = chem.MolFromSmiles(smiles)
        new_order = chem.rdmolfiles.CanonicalRankAtoms(mol)  # reorder the atoms so that they're always in the same canonical order.
        mol = chem.rdmolops.RenumberAtoms(mol, new_order)

        atoms = [atom for atom in mol.GetAtoms()]  # convert PyObject to Julia array
        atoms = sort!(atoms, by = x -> x.GetIdx()+1)  # sort by index to ensure same order as rdkit
        n_atom_feat = sum([length(atom_feat_bins[feat]) for feat in atom_feature_list])
        n_atoms = mol.GetNumAtoms()
        atom_feats = zeros(Float32,n_atoms,n_atom_feat)
        for (a_i,atom) in enumerate(atoms)
            atom_feats[a_i,:] = concat_smiles_atom_feat(smiles_atom_features(atom,feature_list=atom_feature_list))
        end

        bonds = [bond for bond in mol.GetBonds()]  # convert PyObject to Julia array
        bonds = sort!(bonds,by = x->(x.GetBeginAtomIdx()+1,x.GetEndAtomIdx()+1))
        bond_feats = OrderedDict()
        adj_list = [[] for n=1:n_atoms]
        for bond in bonds
            # sort index
            begin_atom_idx = bond.GetBeginAtomIdx()+1
            end_atom_idx = bond.GetEndAtomIdx()+1
            index = begin_atom_idx<end_atom_idx ? (begin_atom_idx,end_atom_idx) : (end_atom_idx, begin_atom_idx)
            bond_feats[index] = concat_smiles_bond_feat(smiles_bond_features(bond,feature_list=bond_feature_list))
            push!(adj_list[begin_atom_idx],end_atom_idx)
            push!(adj_list[end_atom_idx],begin_atom_idx)
        end
        # extract pair features
        bf_len = sum([length(bond_feat_bins[feat]) for feat in bond_feature_list])
        pair_feats = pair_features(mol,bond_feats,adj_list,bf_len = bf_len)
        outputs[i] = WeaveMol(atom_feats, pair_feats, mol.GetNumAtoms(), size(atom_feats,2))
    end
    return outputs
end
