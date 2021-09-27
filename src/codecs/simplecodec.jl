using ..ChemistryFeaturization.AbstractType: AbstractCodec
import ..ChemistryFeaturization.encode
import ..ChemistryFeaturization.decode

export encode, decode

"""
    SimpleCodec(encode_f, decode_f)

SimpleCodec is a most simple, atomic (as in indivisible), and straightforward Codec.

A SimpleCodec is characterized by only the encoding and decoding functions it has as fields.
This means that the encoding-decoding scheme is entirely dependent only upon these functions, and the fields of the FeatureDescriptor to which this Codec is attached. No external parameters (such as number of bins, etc.) are required
by the encoding-decoding implementations.

This Codec is very useful when the encoding-decoding scheme that you want to define for a `FeatureDescriptor` is very simple.
"""
struct SimpleCodec <: AbstractCodec
    encode_f::Function
    decode_f::Function
end

encode(sc::SimpleCodec, val) = sc.encode_f(val)
decode(sc::SimpleCodec, encoded) = sc.decode_f(encoded)