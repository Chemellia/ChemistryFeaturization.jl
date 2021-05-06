# some commentary

# docstring...
mutable struct WeaveMol <: AbstractAtoms
    smiles::String
    elements::Vector{String}
    atom_feats::Vector{SomethingOrOther} # I need to look more carefully to figure this out, heh
    pair_feats::Vector{SomethingOrOther}
    featurization::WeaveFeaturization
    id::Any # probably makes sense to have this in addition to smiles? Maybe?
end

# TODO: add visualize function, either via OpenSMILES or MolecularGraphtion, either via OpenSMILES or MolecularGraph