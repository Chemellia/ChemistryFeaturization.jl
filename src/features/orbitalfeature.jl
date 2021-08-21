include("abstractfeatures.jl")
using DataFrames
using ..ChemistryFeaturization: elements
using ..ChemistryFeaturization.AbstractType: AbstractCodec, AbstractAtoms
using ..ChemistryFeaturization.Codec: SimpleCodec
using ..ChemistryFeaturization.Utils.OrbitalFeatureUtils:
    _name_to_econf, _indexorbital, _econf_to_name
using ..ChemistryFeaturization.Data: valenceshell_conf_df
using SparseArrays

struct OrbitalFeatureDescriptor <: AbstractEnvironmentFeatureDescriptor
    encoder_decoder::AbstractCodec
    function OrbitalFeatureDescriptor()
        # trim down the DataFrame to have only the Symbol, and valence shell configuration (extracted from "Electronic Structure")
        new(SimpleCodec(default_ofd_encode, default_ofd_decode))
    end
end

function (ofd::OrbitalFeatureDescriptor)(a::AbstractAtoms)
    @assert all([el in valenceshell_conf_df[:, :Symbol] for el in elements(a)]) "All elements must be valid and accounted for in the periodic table!"
    ofd.encoder_decoder(ofd, a)
end

# encode
(sc::SimpleCodec)(ofd::OrbitalFeatureDescriptor, a::AbstractAtoms) = sc.encode_f(ofd, a)

# decode
(sc::SimpleCodec)(
    ofd::OrbitalFeatureDescriptor,
    encoded_features::SparseMatrixCSC{Tv,Ti},
) where {Tv,Ti} = sc.decode_f(ofd, encoded_features)

function output_shape(ofd::OrbitalFeatureDescriptor, sc::SimpleCodec)
    # return the number
end

Base.show(io::IO, efd::OrbitalFeatureDescriptor) =
    print(io, "OrbitalFeatureDescriptor created.")

function Base.show(io::IO, ::MIME"text/plain", efd::OrbitalFeatureDescriptor)
    st = "OrbitalFeatureDescriptor created."
    print(io, st)
end

function default_ofd_encode(ofd::OrbitalFeatureDescriptor, a::AbstractAtoms)
    I = Vector{Int16}()
    J = Vector{Int16}()
    V = Vector{Int16}()
    col = 0
    for i in elements(a)
        # X, Y are equivalent to I, V for a SparseVector
        X, Y = _name_to_econf(valenceshell_conf_df, i)
        col += 1    # column number

        for i = 1:length(X)
            append!(I, X[i])
            append!(J, col)
            append!(V, Y[i])
        end
    end
    sparse(I, J, V) # create a SparseMatrix with the valence shell configuration of each element as a column vector
end

function default_ofd_decode(
    ofd::OrbitalFeatureDescriptor,
    encoded_features::SparseArrays.AbstractSparseMatrixCSC{Tv,Ti},
) where {Tv,Ti}
    elements = String[]
    for i = 1:size(encoded_features)[2]
        push!(elements, _econf_to_name(valenceshell_conf_df, encoded_features[:, i]))
    end
    return elements
end
