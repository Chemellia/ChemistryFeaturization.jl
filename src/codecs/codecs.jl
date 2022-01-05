"""
    AbstractCodec

All [Codecs](@ref codecs) defined for different encoding-decoding schemes
must be a subtype of AbstractCodec.
"""
abstract type AbstractCodec end

"""
    encode(codec, val)
Encode `val` according to the scheme described by `codec`.
"""
encode(codec::AbstractCodec, val) = throw(MethodError(encode, codec, val))

"""
    decode(codec, encoded)
Decode `encoded` presuming it was encoded by `codec`.
"""
decode(codec::AbstractCodec, val) = throw(MethodError(decode, codec, val))

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
