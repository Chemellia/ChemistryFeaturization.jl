using ..ChemistryFeaturization.AbstractType: AbstractCodec

"""
    OneHotOneCold(encode_f, decode_f, nbins, logspaced)

AbstractCodec type which uses a dummy variable (as defined in statistical literature), i.e., which employs
one-hot encoding and a one-cold decoding scheme.
"""
struct OneHotOneCold <: AbstractCodec
    encode_f::Function
    decode_f::Function
    nbins::Integer
    logspaced::Bool
end
