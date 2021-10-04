using ..ChemistryFeaturization.AbstractType: AbstractCodec
import ..ChemistryFeaturization.encode
import ..ChemistryFeaturization.decode

"""
    DirectCodec

AbstractCodec type whose encoding function is the identity, or at most some constant * identity.
"""
struct DirectCodec <: AbstractCodec
    scale::Number
end

encode(dc::DirectCodec, val) = dc.scale .* val
decode(dc::DirectCodec, encoded) = encoded ./ dc.scale