"""
    AbstractCodec

All codecs defined for different encoding-decoding schemes
must be a subtype of AbstractCodec.
"""
abstract type AbstractCodec end

"""
    encode(val, codec)
Encode `val` according to the scheme described by `codec`.
"""
encode(val, codec::AbstractCodec) = throw(MethodError(encode, val, codec))

"""
    decode(encoded, codec)
Decode `encoded` presuming it was encoded by `codec`.
"""
decode(val, codec::AbstractCodec) = throw(MethodError(decode, val, codec))

"""
    output_shape(codec)
    output_shape(codec, val)
Return the shape of the encoded output of `codec` (when applied to intput `val`).
"""
output_shape(codec::AbstractCodec) = throw(MethodError(output_shape, codec))
output_shape(codec::AbstractCodec, val) = throw(MethodError(output_shape, codec, val))

include("simplecodec.jl")

include("onehotonecold.jl")

include("directcodec.jl")
