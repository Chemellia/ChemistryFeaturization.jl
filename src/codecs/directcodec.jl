"""
    DirectCodec

Codec type whose encoding function is some constant * the identity function.
"""
struct DirectCodec <: AbstractCodec
    scale::Number
end

encode(dc::DirectCodec, val) = dc.scale .* val
decode(dc::DirectCodec, encoded) = encoded ./ dc.scale
