include("abstractfeatures.jl")

struct PairFeatureDescriptor <: AbstractPairFeatureDescriptor
    name::String
    encode_f::Any
    decode_f::Any
    length::Int # maybe, maybe not (does constrain/assume vector encoding)
    # probably needs some other stuff...
end
