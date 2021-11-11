include("abstractfeatures.jl")

"""
    BondFeatureDescriptor{A}

A feature object that encodes features associated with bonded pairs of atoms in a structure. For pair features, see `PairFeatureDescriptor`.

Type parameter represents the structure representation(s) from which this feature descriptor is able to compute features.

## Fields
- `name::String`: the name of the feature
- `compute_f`: function that takes in a structure of type `A` and returns values of this feature for every bonded pair of atoms in the structure
- `encoder_decoder::AbstractCodec`: codec that encodes/decodes values of this feature
- `categorical::Bool`: flag for whether the feature is categorical or continuous-valued
- `encodable_elements::Vector{String}`: list of elements (by symbol) that can be encoded by this feature
"""
struct BondFeatureDescriptor{A,C<:AbstractCodec} <: AbstractPairFeatureDescriptor
    name::String
    compute_f::Function
    encoder_decoder::C
    categorical::Bool
    encodable_elements::Vector{String}
end

function get_value(bfd::BondFeatureDescriptor{GraphMol}, a::AbstractAtoms{GraphMol})
    @assert all([el in encodable_elements(bfd) for el in elements(a)]) "Feature $(bfd.name) cannot encode some element(s) in this structure!"
    vals = bfd.compute_f(a.structure)
    bonds = a.structure.edges
    mat = Matrix{Union{Missing,eltype(vals)}}(missing, length(a.structure.nodeattrs), length(a.structure.nodeattrs))
    for i in 1:length(bonds)
        mat[bonds[i][1], bonds[i][2]] = vals[i]
        mat[bonds[i][2], bonds[i][1]] = vals[i]
    end
    return mat
end
