module ChemistryFeaturization

export inverse_square, exp_decay, build_graph, visualize_graph
include("build_cif_graphs.jl")
include("graph_vis.jl")

export atom_data_df, make_feature_vectors, decode_feature_vector
include("species_features.jl")

# Sean: add which functions/variables to be exported from your files; I've just added the ones to make the tests pass
export weave_featurize, smiles_atom_features, smiles_bond_features
include("smiles_rdkit_features.jl")

end