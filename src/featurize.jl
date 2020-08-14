# this file will have the "master" function that's actually exposed to the user
include("cif_fcns.jl")
include("smiles_fcns.jl")

function featurize_inputs(input_data_folder, SMILES_atom_features, SMILES_bond_features, other_atom_features, graphs_folder, features_folder)
    # determine type of data in input folder (must be all SMILES or all CIFs)
    # check whether user has asked for features that make sense (may require some ancillary fcns)
    # loop over input files and farm out actual parsing to specialized functions for CIF vs. SMILES
    # eventually, save out adjacency matrices to graphs_folder and feature vectors/metadata to features_folder
end
