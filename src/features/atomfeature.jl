using ..ChemistryFeaturization.Utils.AtomFeatureUtils
using DataFrames

#=
Feature of a single atom.

May be contextual (depends on neighborhood) or elemental (defined just by the atomic identity of the node).
=#
# proper docstring
# TODO: add way to get range/list of possible values for feature...
struct AtomFeature <: AbstractFeature
    name::String
    encode_f::Any
    decode_f::Any
    categorical::Bool
    contextual::Bool # can get from elemental lookup table (false) or not (true)?
    length::Int # length of encoded vector
    encodable_elements::Vector{String}
end

# pretty printing, short version
Base.show(io::IO, af::AtomFeature) = print(io, "AtomFeature $(af.name)")

# pretty printing, long version
function Base.show(io::IO, ::MIME"text/plain", af::AtomFeature)
    st = "AtomFeature $(af.name):\n   categorical: $(af.categorical)\n   contextual: $(af.contextual)\n   encoded length: $(af.length)"
    print(io, st)
end

# docstring
function AtomFeature(
    feature_name;
    nbins = default_nbins,
    lookup_table::DataFrame = atom_data_df,
    logspaced::Bool = default_log(feature_name, lookup_table),
    categorical::Bool = default_categorical(feature_name, lookup_table),
)
    colnames = names(lookup_table)
    @assert feature_name in colnames && "Symbol" in colnames "Your lookup table must have a column called :Symbol and one with the same name as your feature to be usable!"

    local vector_length
    if categorical
        vector_length = length(unique(skipmissing(lookup_table[:, Symbol(feature_name)])))
    else
        vector_length = nbins
    end
    encode_f =
        atoms -> reduce(
            hcat,
            map(
                e -> onehot_lookup_encoder(
                    e,
                    feature_name;
                    lookup_table = lookup_table,
                    nbins = nbins,
                    logspaced = logspaced,
                    categorical = categorical,
                    
                ),
                atoms.elements,
            ),
        )
    decode_f =
        encoded -> onecold_decoder(
            encoded,
            feature_name;
            lookup_table = lookup_table,
            nbins = nbins,
            logspaced = logspaced,
            categorical = categorical,
        )
    AtomFeature(
        feature_name,
        encode_f,
        decode_f,
        categorical,
        false,
        vector_length,
        ChemistryFeaturization.Utils.AtomFeatureUtils.encodable_elements(
            feature_name,
            lookup_table,
        ),
    )
end

encodable_elements(f::AtomFeature) = f.encodable_elements

function (f::AtomFeature)(a::AbstractAtoms)
    @assert all([el in f.encodable_elements for el in a.elements]) "Feature $(f.name) cannot encode some element(s) in this structure!"
    f.encode_f(a)
end

# TODO: some Weave stuff needed here
