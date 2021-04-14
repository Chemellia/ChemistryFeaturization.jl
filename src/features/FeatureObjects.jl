#= FEATURE OBJECTS

I decided on an abstract type to allow the dispatch shown below in the "magical encoding"
and "magical decoding" bits.

Tn is the "natural" type of feature values, Te is the type of encoded values for the WHOLE
STRUCTURE (not e.g. for a single atom).

Representative examples – feature: feature object type, Tn, Te:

block (s,p,d,f):            AtomFeat (contextual=false), String, Flux.OneHotMatrix 
electronegativity (binned): AtomFeat (contextual=false), Float32, Flux.OneHotMatrix
electronegativity (direct): AtomFeat (contextual=false), Float32, Vector{Float32}
oxidation state:            AtomFeat (contextual=true) Int, Vector{Float32}
distance between atoms:     PairFeat, Float32, Matrix{Float32}
bond type:                  PairFeat, String, Array{Float32,3}

(bond type is categorical and gets one-hot encoded, so the first two indices of the Matrix
should be atom indices and the third should be indexing into the one-hot vector, which should just contain all zeros if the two atoms are not bonded (or we add a bin to the one-hot encoding to indicate that)

All subtypes should define `encode_f` and `decode_f`
`encode_f` should take in an atoms object and return Te
`decode_f` should take in something of type Te and return something of type Tn
=#

#module FeatureObjects

# link to guidance in docs about how to implement new feature types

#using .ChemistryFeaturization:

# export...
export AtomFeat, PairFeat

# include...
include("atomfeat.jl")
include("pairfeat.jl")

# generic encode
# docstring
function (f::AbstractFeature{Tn,Te})(a::AbstractAtoms) where {Te,Tn}
    f.encode_f(a)
end

# generic decode
# docstring
function decode(f::AbstractFeature{Tn,Te}, encoded_f::Te) where {Te,Tn}
    f.decode_f(encoded_f)
end

#end
