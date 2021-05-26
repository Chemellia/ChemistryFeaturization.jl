using ..ChemistryFeaturization.Utils.AtomFeatureUtils
using DataFrames

# TODO: figure out what scheme would look like that is flexible to direct-value encoding (may just need a different feature type since it'll have to handle normalization, etc. too)
"""
    ElementFeature(feature_name, encode_f, decode_f, categorical, contextual, length, encodable_elements)

Construct a feature object that encodes features associated with individual atoms that depend only upon their elemental identity (if you want to encode a feature that depends upon an atom's environment, you shold use SpeciesFeature!)

## Arguments
- `name::String`: the name of the feature
- `categorical::Bool`: flag for whether the feature is categorical or continuous-valued
- `length::Int`: length of encoded vector
- `logspaced::Bool`: whether onehot-style bins should be logarithmically spaced or not
- `lookup_table::DataFrame`: table containing values of feature for every encodable element
"""
struct ElementFeature <: AbstractAtomFeature
    name::String
    categorical::Bool
    length::Integer
    logspaced::Bool
    lookup_table::DataFrame
end

# TODO: update this, encoder stuff needs to be broken out as dispatches
# also, should trim lookup_table to just have the columns it needs before constructing object
function ElementFeature(
    feature_name::String,
    lookup_table::DataFrame = atom_data_df;
    length::Integer = default_nbins,
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
                    feature_name,
                    lookup_table;
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
            feature_name,
            lookup_table;
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
        encodable_elements(feature_name, lookup_table),
    )
end

# TODO: update this to just call thing below essentially
#encodable_elements(f::ElementFeature) = f.encodable_elements

function encodable_elements(feature_name::String, lookup_table::DataFrame = atom_data_df)
    info = lookup_table[:, [Symbol(feature_name), :Symbol]]
    return info[
        findall(x -> !ismissing(x), getproperty(info, Symbol(feature_name))),
        :Symbol,
    ]
end