#=
Feature of a pair of atoms. Currently only used in WeaveModel, but may eventually
have a version of AtomGraph that allows edge features...
=#
struct PairFeature<:AbstractFeature
    name
    encode_f
    decode_f
    length::Int # maybe, maybe not (does constrain/assume vector Te)
    # probably needs some other stuff...
end

# copy/translate weave stuff
# a note: I think it's more elegant (and probably efficient?) if pair features assert
# always arrays where the first two indices are atom indices and any further are indexing
# into the feature if it is itself higher-dimensional. This is not how the weave featurizer
# currently works...