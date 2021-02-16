module ChemistryFeaturization

export AtomFeat, atom_data_df, build_atom_feats, make_feature_vectors, decode_feature_vector, default_nbins
include("atomfeat.jl")

export AtomGraph, normalized_laplacian, add_features!, add_features_batch!, visualize_graph
include("atomgraph.jl")

export inverse_square, exp_decay, build_graph, build_graphs_batch, read_graphs_batch
include("pmg_graphs.jl")
using .graph_building: inverse_square, exp_decay, build_graph, build_graphs_batch, read_graphs_batch

# TODO: possibly move all the Weave stuff to another package altogether, if not need to tidy up modules/exports
export smiles_atom_features, smiles_bond_features
include("weave_fcns.jl")
using .weave_fcns: smiles_atom_features, smiles_bond_features, chem

export weave_featurize
include("featurize.jl")
using .featurize: weave_featurize

end
