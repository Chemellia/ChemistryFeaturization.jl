#=
Feature of a single atom.

May be contextual (depends on neighborhood) or elemental (defined just by the atomic identity of the node).
=#
# proper docstring
struct AtomFeature<:AbstractFeature
    name::String
    encode_f
    decode_f
    categorical::Bool
    contextual::Bool # can get from elemental lookup table (false) or not (true)?
    length::Int # length of encoded vector
end

# TODO: add pretty printing

# NEXT: test this constructor
# docstring
using ..ChemistryFeaturization.Utils.AtomFeatureUtils
function AtomFeature(feature_name; nbins=default_nbins, logspaced=false)
    @assert feature_name in continuous_feature_names || feature_name in categorical_feature_names "Cannot automatically build AtomFeat for $feature_name; I can't find it in a lookup table!"
    local vector_length
    categorical = feature_name in categorical_feature_names
    if categorical
        vector_length = length(categorical_feature_vals[feature_name])
    else
        vector_length = nbins
    end
    encode_f = atoms -> map(e->onehot_lookup_encoder(e, feature_name; nbins=nbins, logspaced=logspaced), atoms.elements)
    decode_f = encoded -> onecold_lookup_decoder(encoded, feature_name; nbins=nbins, logspaced=logspaced)
    AtomFeature(feature_name, encode_f, decode_f, categorical, false, vector_length)
end

#=
we'll define a bunch of automatic stuff for building AtomFeatures with built-in data
that will essentially copy the current contents of atomfeat.jl, ideally with some
additions enabling the user to augment the lookup table

should also have convenience functions for building encode_f and decode_f via keywords for, e.g.:
    choosing one-hot (and how many bins, etc.) vs. direct float encoding
    logspaced vs. linear spaced binning for one-hot encoding
    maybe other stuff too
...
=#

# then basically copy a bunch of the weave stuff, plus potentially things like oxidation state etc. from pymatgen belong here too