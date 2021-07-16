include("abstractfeatures.jl")
using DataFrames
using ..ChemistryFeaturization: elements
using ..ChemistryFeaturization.AbstractType: AbstractCodec, AbstractAtoms
using ..ChemistryFeaturization.Codec: SimpleCodec
using ..ChemistryFeaturization.Utils.OrbitalFeatureUtils:
    valence_shell_config, _indexorbital
using SparseArrays

struct OrbitalFeatureDescriptor <: AbstractEnvironmentFeatureDescriptor
    lookup_table::DataFrame
    encoder_decoder::AbstractCodec
    OrbitalFeatureDescriptor() =
        new(atom_data_df, SimpleCodec(default_ofd_encode, default_ofd_decode))
end

function (ofd::OrbitalFeatureDescriptor)(a::AbstractAtoms)
    @assert all([el in ofd.lookup_table[:, :Symbol] for el in elements(a)]) "All elements must be valid and accounted for in the periodic table!"
    ofd.encoder_decoder(ofd, a)
end

# encode
(sc::SimpleCodec)(ofd::OrbitalFeatureDescriptor, a::AbstractAtoms) = sc.encode_f(ofd, a)

# decode
(sc::SimpleCodec)(ofd::OrbitalFeatureDescriptor, encoded_features::SparseMatrixCSC{Tv,Ti},
) where {Tv,Ti} = sc.decode_f(encoded_features)

function output_shape(ofd::OrbitalFeatureDescriptor, sc::SimpleCodec)
    # return the number
end

function default_ofd_encode(ofd::OrbitalFeatureDescriptor, a::AbstractAtoms)
    I = Vector{Int16}()
    J = Vector{Int16}()
    V = Vector{Int16}()
    col = 0
    for i in elements(a)
        X, Y = valence_shell_config(ofd.lookup_table, i)
        col += 1

        for i = 1:length(X)
            append!(I, X[i])
            append!(J, col)
            append!(V, Y[i])
        end
    end
    sparse(I, J, V)
end

function default_ofd_decode(
    encoded_features::SparseArrays.AbstractSparseMatrixCSC{Tv,Ti},
) where {Tv,Ti}
    elements = String[]
    rows = rowvals(encoded_features)
    vals = nonzeros(encoded_features)
    _, n = size(encoded_features)
    shell_conf = ""
    for j = 1:n
        for i in nzrange(encoded_features, j)
            row = rows[i]
            val = vals[i]
            shell_conf = shell_conf * _indexorbital(row) * string(val) * "."
        end
        element = (filter(
            "Valence Shell Configuration" => x -> x == chop(shell_conf),
            atom_data_df;
            view = true,
        )[
            !,
            :Symbol,
        ])[1]
        append!(elements, (element,))
        shell_conf = ""
    end
    return elements
end
