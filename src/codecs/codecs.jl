"""
    Codec

Codecs provide the encoding-decoding scheme in a flexible way.
Codecs are titled according to the general attributes that the codec scheme
would require, and NOT the actual encoding/decoding function itself.

This makes things customizable, and allows plug-and-play behaviour with
different variants of the same codec scheme, or for strikingly similar codec
schemes.
"""
module Codec

include("simplecodec.jl")
export SimpleCodec

include("onehotonecold.jl")
export OneHotOneCold

export encode, decode

end
