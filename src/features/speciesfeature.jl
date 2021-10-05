using DataFrames
using ..ChemistryFeaturization.Utils.SpeciesFeatureUtils
using ..ChemistryFeaturization.Codec

include("abstractfeatures.jl")

"""
    SpeciesFeatureDescriptor{A}

A feature object that encodes features associated with individual atoms that depend upon their local environment in some way (if your feature is defined only by elemental identity, you should use ElementFeatureDescriptor!)

Type parameter represents the structure representation(s) from which this feature descriptor is able to compute features.

## Fields
- `name::String`: the name of the feature
- `compute_f`: function that takes in a structure of type `A` and returns values of this feature for every atom in the structure
- `encoder_decoder::AbstractCodec`: codec that encodes/decodes values of this feature
- `categorical::Bool`: flag for whether the feature is categorical or continuous-valued
- `encodable_elements::Vector{String}`: list of elements (by symbol) that can be encoded by this feature
"""
struct SpeciesFeatureDescriptor{A} <: AbstractAtomFeatureDescriptor
    name::String
    compute_f::Any
    encoder_decoder::AbstractCodec
    categorical::Bool
    encodable_elements::Vector{String}
end

"""
    SpeciesFeatureDescriptor(name::String)
"""
function SpeciesFeatureDescriptor(name::String)
    @assert name in keys(sfd_names_props) "That species feature isn't one of the built-in ones, you'll have to construct it directly. Consult ChemistryFeaturiztion.Utils.SpeciesFeatureUtils.sfd_names_props for built-in species features."

    info = sfd_names_props[name]
    categorical = info[:categorical]
    local codec
    if categorical
        possible_vals = info[:possible_vals]
        if length(possible_vals) == 2
            @assert possible_vals == [true, false]
            codec = DirectCodec(1)
        else
            codec = OneHotOneCold(true, possible_vals)
        end
    else
        # TODO: figure out default binning situation for continuous-valued SFD's
        #codec = OneHotOneCold(false, )
    end
    SpeciesFeatureDescriptor{info[:A]}(
        name,
        info[:compute_f],
        codec,
        categorical,
        info[:encodable_elements],
    )
end

# pretty printing, short version
Base.show(io::IO, af::SpeciesFeatureDescriptor) = print(io, "SpeciesFeature $(af.name)")

# pretty printing, long version
function Base.show(io::IO, ::MIME"text/plain", af::SpeciesFeatureDescriptor{A}) where {A}
    st = "SpeciesFeature $(af.name):\n   categorical: $(af.categorical)\n   works on: $(A)"
    print(io, st)
end

# TODO: add way to get range/list of possible values for feature...
encodable_elements(f::SpeciesFeatureDescriptor) = f.encodable_elements


function get_value(sfd::SpeciesFeatureDescriptor{A}, a::AbstractAtoms{<:A}) where {A}
    @assert all([el in encodable_elements(sfd) for el in elements(a)]) "Feature $(efd.name) cannot encode some element(s) in this structure!"
    sfd.compute_f(a.structure)
end
