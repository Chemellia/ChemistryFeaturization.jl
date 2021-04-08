# some commentary

# docstring...
mutable struct WeaveMol <: AbstractAtoms
    smiles::String
    elements::Vector{String}
    features::Tuple{SomethingOrOther} # I need to look more carefully to figure this out heh
    featurization::WeaveFeaturization
    id # probably makes sense to have this in addition to smiles? Maybe?
end