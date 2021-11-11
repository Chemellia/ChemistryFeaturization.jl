include("abstractfeatures.jl")

"""
    PairFeatureDescriptor{A}

A feature object that encodes features associated with pairs of atoms in a structure that may or may not be bonded to each other. For bond features, see `BondFeatureDescriptor`.

Type parameter represents the structure representation(s) from which this feature descriptor is able to compute features.

## Fields
- `name::String`: the name of the feature
- `compute_f`: function that takes in a structure of type `A` and returns values of this feature for every pair of atoms in the structure organized in a matrix
- `encoder_decoder::AbstractCodec`: codec that encodes/decodes values of this feature
- `categorical::Bool`: flag for whether the feature is categorical or continuous-valued
- `encodable_elements::Vector{String}`: list of elements (by symbol) that can be encoded by this feature
"""
struct PairFeatureDescriptor{A,C<:AbstractCodec} <: AbstractPairFeatureDescriptor
    name::String
    compute_f::Function
    encoder_decoder::C
    categorical::Bool
    encodable_elements::Vector{String}
end

function get_value(pfd::PairFeatureDescriptor{A}, a::AbstractAtoms{<:A}) where {A}
    @assert all([el in encodable_elements(pfd) for el in elements(a)]) "Feature $(pfd.name) cannot encode some element(s) in this structure!"
    pfd.compute_f(a.structure)
end
