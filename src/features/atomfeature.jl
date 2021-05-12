using ..ChemistryFeaturization.Utils.AtomFeatureUtils

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
    logspaced = feature_name in keys(default_log) ? default_log[feature_name] : false,
    categorical = feature_name in categorical_feature_names,
    lookup_table = atom_data_df,
)
    local builtin =
        feature_name in continuous_feature_names ||
        feature_name in categorical_feature_names
    if !builtin
        @assert feature_name in names(lookup_table) "$feature_name is not a built-in feature, but you haven't provided a lookup table to find its values!"
    end
    local vector_length
    if categorical
        if builtin
            vector_length = length(categorical_feature_vals[feature_name])
        else
            vector_length = length(unique(lookup_table[:, Symbol(feature_name)]))
        end
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
                    nbins = nbins,
                    logspaced = logspaced,
                    categorical = categorical,
                    lookup_table = lookup_table,
                ),
                atoms.elements,
            ),
        )
    decode_f =
        encoded -> onecold_decoder(
            encoded,
            feature_name;
            nbins = nbins,
            logspaced = logspaced,
            categorical = categorical,
            lookup_table = lookup_table,
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
