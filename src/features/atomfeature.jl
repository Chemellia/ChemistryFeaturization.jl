#=
Feature of a single atom.

May be contextual (depends on neighborhood) or elemental (defined just by the atomic identity of the node).
=#
# proper docstring
struct AtomFeature{Tn,Te}<:AbstractFeature{Tn,Te}
    name::String
    encode_f
    decode_f
    categorical::Bool
    contextual::Bool # can get from elemental lookup table (false) or not (true)?
    length::Int # length of encoded vector
end

# TODO: add pretty printing

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