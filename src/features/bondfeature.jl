using ..ChemistryFeaturization.Utils.BondFeatureUtils
using ..ChemistryFeaturization.Codec
using ..ChemistryFeaturization.Atoms

using MolecularGraph

include("abstractfeatures.jl")

"""
    BondFeatureDescriptor{A}

A feature object that encodes features associated with bonded pairs of atoms in a structure. For general pair features, see `PairFeatureDescriptor`.

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

# pretty printing, long version
function Base.show(io::IO, ::MIME"text/plain", fd::BondFeatureDescriptor{A}) where {A}
    st = "$(typeof(fd)) $(fd.name):\n   categorical: $(fd.categorical)\n   works on: $(A)"
    print(io, st)
end

"""
    BondFeatureDescriptor(name::String)

Construct a `BondFeatureDescriptor` from one of the built-in options (currently supports a variety of MolecularGraph functions, see Utils.BondFeatureUtils.bfd_names_props for details).
"""
function BondFeatureDescriptor(name::String)
    @assert name in keys(bfd_names_props) "That species feature isn't one of the built-in ones, you'll have to construct it directly. Consult ChemistryFeaturiztion.Utils.BondFeatureUtils.bfd_names_props for built-in species features."

    info = bfd_names_props[name]
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
        # TODO: figure out default binning situation for continuous-valued BFD's
        #codec = OneHotOneCold(false, )
    end
    BondFeatureDescriptor{info[:A],typeof(codec)}(
        name,
        info[:compute_f],
        codec,
        categorical,
        info[:encodable_elements],
    )
end
# for some reason having AbstractAtoms instead of AtomGraph doesn't get dispatched propertly, I assume this has something to do with the type parameters...
function get_value(
    bfd::BondFeatureDescriptor{GraphMol},
    a::AtomGraph{GraphMol{A,B}},
) where {A,B}
    @assert all([el in encodable_elements(bfd) for el in elements(a)]) "Feature $(bfd.name) cannot encode some element(s) in this structure!"
    vals = bfd.compute_f(a.structure)
    bonds = a.structure.edges
    mat = Matrix{Union{Missing,eltype(vals)}}(
        missing,
        length(a.structure.nodeattrs),
        length(a.structure.nodeattrs),
    )
    # force it to be symmetric...TODO: add option for asymmetric features
    for i = 1:length(bonds)
        mat[bonds[i][1], bonds[i][2]] = vals[i]
        mat[bonds[i][2], bonds[i][1]] = vals[i]
    end
    return mat
end
