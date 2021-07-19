"""
    Codec

Codecs provide the encoding-decoding scheme in a flexible way.
Codecs are titled according to the general attributes that the codec scheme
would require, and NOT the actual encoding/decoding function itself.

Every Codec object MUST have one set of `encode_f::Function` and
`decode_f::Function` fields.

This makes things customizable, and allows plug-and-play behaviour with
different variants of the same codec scheme, or for strikingly similar codec
schemes.
"""
module Codec

@enum EncodeOrDecode ENCODE DECODE
export EncodeOrDecode

include("OneHotOneCold.jl")
export OneHotOneCold, build_onehot_vec

end
