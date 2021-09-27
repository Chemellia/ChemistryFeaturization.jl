include("abstractfeatures.jl")
using DataFrames
using ..ChemistryFeaturization: elements
using ..ChemistryFeaturization.AbstractType: AbstractCodec, AbstractAtoms
using ..ChemistryFeaturization.Codec: SimpleCodec
using ..ChemistryFeaturization.Utils.OrbitalFeatureUtils:
    _indexorbital, _econf_to_name, _orbitalsparse
using ..ChemistryFeaturization.Data: valenceshell_conf_df
using SparseArrays

"""
    OrbitalFeatureDescriptor

Describes the orbital configuration of an element.
"""
struct OrbitalFeatureDescriptor <: AbstractEnvironmentFeatureDescriptor
    encoder_decoder::AbstractCodec
    function OrbitalFeatureDescriptor()
        # trim down the DataFrame to have only the Symbol, and valence shell configuration (extracted from "Electronic Structure")
        new(SimpleCodec(default_ofd_encode, default_ofd_decode))
    end
end

function get_value(ofd::OrbitalFeatureDescriptor, el::String)
    @assert el in valenceshell_conf_df[:, :Symbol] "All elements must be valid and accounted for in the periodic table!"

    getproperty(valenceshell_conf_df[valenceshell_conf_df.Symbol.==el, :][1,:], Symbol("Electronic Structure"))
end

get_value(ofd::OrbitalFeatureDescriptor, a::AbstractAtoms) = map(e -> get_value(ofd, e), elements(a))

(ofd::OrbitalFeatureDescriptor)(el::String) = get_value(ofd, el)

function output_shape(ofd::OrbitalFeatureDescriptor, sc::SimpleCodec)
    # return the number
end

Base.show(io::IO, efd::OrbitalFeatureDescriptor) =
    print(io, "OrbitalFeatureDescriptor created.")

function Base.show(io::IO, ::MIME"text/plain", efd::OrbitalFeatureDescriptor)
    st = "OrbitalFeatureDescriptor created."
    print(io, st)
end

"""
    default_ofd_encode(ofd::OrbitalFeatureDescriptor, a::AbstractAtoms)

Default encoding scheme for OrbitalFeatureDescriptor, which returns a
sparse matrix, where column_i is a sparse representation of the electronic
configuration of the ith element in the [Atoms](@ref atoms) object.
"""
function default_ofd_encode(elec_configs::Vector{String})
    I = Vector{Int16}()
    J = Vector{Int16}()
    V = Vector{Int16}()
    col = 0
    for s in elec_configs
        # X, Y are equivalent to I, V for a SparseVector
        X, Y = _orbitalsparse(s)
        col += 1    # column number

        for i = 1:length(X)
            append!(I, X[i])
            append!(J, col)
            append!(V, Y[i])
        end
    end
    sparse(I, J, V) # create a SparseMatrix with the valence shell configuration of each element as a column vector
end


"""
    default_ofd_decode(
        ofd::OrbitalFeatureDescriptor,
        encoded_features::SparseArrays.AbstractSparseMatrixCSC{Tv,Ti}
    ) where {Tv,Ti}

Default decoding scheme for OrbitalFeatureDescriptor, which returns the 
vector of elements encoded, given the encoded sparse matrix.
"""
function default_ofd_decode(encoded_features::SparseArrays.AbstractSparseMatrixCSC{Tv,Ti}) where {Tv,Ti}
    elements = String[]
    for i = 1:size(encoded_features)[2]
        push!(elements, _econf_to_name(encoded_features[:, i]))
    end
    return elements
end
