using DataFrames

"""
    SpeciesFeatureDescriptor(feature_name, encode_f, decode_f, categorical, contextual, length, encodable_elements)

Construct a feature object that encodes features associated with individual atoms that depend upon their local environment in some way (if your feature is defined only by elemental identity, you should use ElementFeatureDescriptor!)

## Arguments
- `name::String`: the name of the feature
- `encode_f::Function`: a function that takes in <:AbstractAtoms and returns encoded values of this feature for the atoms in that structure
- `decode_f::Function`: inverse function to `encode_f`, takes in encoded feature and returns value (for categorical) or range of values (for continuous-valued) of the feature
- `categorical::Bool`: flag for whether the feature is categorical or continuous-valued
- `length::Int`: length of encoded vector
- `encodable_elements::Vector{String}`: list of elements (by symbol) that can be encoded by this feature
"""
struct SpeciesFeatureDescriptor <: AbstractAtomFeatureDescriptor
    name::String
    encode_f::Function
    decode_f::Function
    categorical::Bool
    length::Integer
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

function (f::SpeciesFeatureDescriptor)(a::AbstractAtoms)
    @assert all([el in f.encodable_elements for el in a.elements]) "Feature $(f.name) cannot encode some element(s) in this structure!"
    f.encode_f(a)
end

# TODO: some Weave stuff needed here?
