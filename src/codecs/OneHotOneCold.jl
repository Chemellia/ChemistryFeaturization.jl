using ..ChemistryFeaturization.AbstractTypes: Codec

"""
    OneHotOneCold(encode_f, decode_f, nbins, logspaced)

Codec type which uses a dummy variable (as defined in statistical literature), i.e., which employs
one-hot encoding and a one-cold decoding scheme.
"""
struct OneHotOneCold <: Codec
    encode_f::Function
    decode_f::Function
    nbins::Integer
    logspaced::Bool
end

