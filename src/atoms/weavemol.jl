# some commentary

# docstring...

mutable struct WeaveMol <: AbstractAtoms # maybe we want to add AbstractAtomGraph <: AbstractAtoms that this and AtomGraph would subtype?
    smiles::String
    elements::Vector{String}
    graph::SimpleWeightedGraph{<:Integer,<:Real}
    encoded_atom_features::Union{Matrix{<:Real},Nothing}
    encoded_pair_features::Union{Matrix{Union{Vector{<:Real},Nothing}},Nothing}
    featurization::WeaveFeaturization
    id::Any # probably makes sense to have this in addition to smiles? Maybe?
end
