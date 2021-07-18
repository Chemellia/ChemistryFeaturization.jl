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
        # X, Y are equivalent to I, V for a SparseVector
        X, Y = valence_shell_config(ofd.lookup_table, i)
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
    encoded_features::SparseArrays.AbstractSparseMatrixCSC{Tv,Ti},
) where {Tv,Ti}
    elements = String[]
    rows = rowvals(encoded_features)  # the rows of the sparse matrix
    vals = nonzeros(encoded_features) # the values in the sparse matrix
    m, n = size(encoded_features) # 'dimensions' of the sparse matrix
    # for each column vector (remember, a column vector is representative of an element's valence shell configuration)
    for j = 1:n 
        shell_conf = "" # the shell configuration for an element
        for i in nzrange(encoded_features, j)   # for each valid value in a column
            row = rows[i]  
            val = vals[i]
            # find the configuration for each shell and concatenate it to the `shell_conf`
            shell_conf = shell_conf * _indexorbital(row) * string(val) * "." 
        end
        # `chop` the trailing '.', and rearrange order of the shells to match the format followed in the default `lookup_table`
        shell_conf = join(sort(split(chop(shell_conf), '.'), by=first), '.')

        # find the DataFrame that has the matching `shell_conf` as its Valence Shell Configuration, and get its `Symbol`
        element = (filter(
            "Valence Shell Configuration" => x -> x == shell_conf,
            atom_data_df;
            view = true,
        )[
            !,
            :Symbol,
        ])[1]
        # append the `Symbol` that was just extracted to the container meant for storing all the elements in the encoded Atoms object
        append!(elements, (element,))
    end
    return elements
end
