# some commentary

# docstring...

# represents a given pair of atoms by their respective indices in the adjacency matrix representation
const PairIndexes = Tuple{Integer,Integer}

mutable struct WeaveMol <: AbstractAtoms
    smiles::String
    elements::Vector{String}
    atom_features::Vector{SomethingOrOther} # Details TBD
    # pairs encapsulates the WeavePair features for a given pair of atoms in a molecule represented using PairIndexes
    pairs::Vector{Tuple{PairIndexes,WeavePair}}
    featurization::WeaveFeaturization
    # graph::SimpleWeightedGraph{<:Integer,<:Real} # is something like an adjacency graph required here?
    id::Any # probably makes sense to have this in addition to smiles? Maybe?
end

# TODO: add visualize function, either via OpenSMILES or MolecularGraphtion, either via OpenSMILES or MolecularGraph
