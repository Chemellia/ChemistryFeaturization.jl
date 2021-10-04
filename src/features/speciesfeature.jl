using DataFrames

include("abstractfeatures.jl")

"""
    SpeciesFeatureDescriptor

Construct a feature object that encodes features associated with individual atoms that depend upon their local environment in some way (if your feature is defined only by elemental identity, you should use ElementFeatureDescriptor!)

Type parameter represents the structure representation(s) from which this feature descriptor is able to compute features.

## Fields
- `name::String`: the name of the feature
- `encode_f::Function`: a function that takes in <:AbstractAtoms and returns encoded values of this feature for the atoms in that structure
- `decode_f::Function`: inverse function to `encode_f`, takes in encoded feature and returns value (for categorical) or range of values (for continuous-valued) of the feature
- `categorical::Bool`: flag for whether the feature is categorical or continuous-valued
- `length::Int`: length of encoded vector
- `encodable_elements::Vector{String}`: list of elements (by symbol) that can be encoded by this feature
"""
struct SpeciesFeatureDescriptor{A} <: AbstractAtomFeatureDescriptor
    name::String
    compute_f::Any
    encoder_decoder::AbstractCodec
    categorical::Bool
    encodable_elements::Vector{String}
end


# pretty printing, short version
Base.show(io::IO, af::SpeciesFeatureDescriptor) = print(io, "AtomFeature $(af.name)")

# pretty printing, long version
function Base.show(io::IO, ::MIME"text/plain", af::SpeciesFeatureDescriptor)
    st = "AtomFeature $(af.name):\n   categorical: $(af.categorical)\n   encoded length: $(af.length)"
    print(io, st)
end

# TODO: add way to get range/list of possible values for feature...
encodable_elements(f::SpeciesFeatureDescriptor) = f.encodable_elements


function get_value(sfd::SpeciesFeatureDescriptor{A}, a::AbstractAtoms{A}) where {A}
    @assert all([el in encodable_elements(efd) for el in elements(a)]) "Feature $(efd.name) cannot encode some element(s) in this structure!"
    sfd.compute_f(a.structure) # TODO: currently inconsistent with analogous behavior of EFD
end
