module ChemistryFeaturization

export AtomFeat, atom_data_df, build_atom_feats, make_feature_vectors, decode_feature_vector, default_nbins
include("atomfeat.jl")

export AtomGraph, normalized_laplacian, add_features!, visualize_graph
include("atomgraph.jl")

export inverse_square, exp_decay, build_graph, build_graphs_batch
include("pmg_graphs.jl")

# Sean: add which functions/variables to be exported from your files; I've just added the ones to make the tests pass
export weave_featurize, smiles_atom_features, smiles_bond_features
include("smiles_rdkit_features.jl")

end