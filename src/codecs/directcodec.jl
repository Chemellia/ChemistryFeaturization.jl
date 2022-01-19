"""
    DirectCodec

Codec type whose encoding function is some constant * the identity function.
"""
struct DirectCodec <: AbstractCodec
    scale::Number
end

encode(val, dc::DirectCodec) = dc.scale .* val
decode(encoded, dc::DirectCodec) = encoded ./ dc.scale
