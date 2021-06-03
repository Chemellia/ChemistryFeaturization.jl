using ..ChemistryFeaturization.Utils.ElementFeatureUtils
using DataFrames

# TODO: figure out what scheme would look like that is flexible to direct-value encoding (may just need a different feature type since it'll have to handle normalization, etc. too)
"""
    ElementFeatureDescriptor(feature_name, encode_f, decode_f, categorical, contextual, length, encodable_elements)

Construct a feature object that encodes features associated with individual atoms that depend only upon their elemental identity (if you want to encode a feature that depends upon an atom's environment, you shold use SpeciesFeatureDescriptor!)

## Arguments
- `name::String`: the name of the feature
- `categorical::Bool`: flag for whether the feature is categorical or continuous-valued
- `length::Int`: length of encoded vector
- `logspaced::Bool`: whether onehot-style bins should be logarithmically spaced or not
- `lookup_table::DataFrame`: table containing values of feature for every encodable element
"""
struct ElementFeatureDescriptor <: AbstractAtomFeatureDescriptor
    name::String
    length::Integer
    nbins::Integer
    logspaced::Bool
    categorical::Bool
    lookup_table::DataFrame
end

# TODO: update this, encoder stuff needs to be broken out as dispatches
# also, should trim lookup_table to just have the columns it needs before constructing object
function ElementFeatureDescriptor(
    feature_name::String,
    lookup_table::DataFrame = atom_data_df;
    nbins::Integer = default_nbins,
    logspaced::Bool = default_log(feature_name, lookup_table),
    categorical::Bool = default_categorical(feature_name, lookup_table)
)
    colnames = names(lookup_table)
    @assert feature_name in colnames && "Symbol" in colnames "Your lookup table must have a column called :Symbol and one with the same name as your feature to be usable!"

    local vector_length
    if categorical
        vector_length = length(unique(skipmissing(lookup_table[:, Symbol(feature_name)])))
    else
        vector_length = nbins
    end

    lookup_table = lookup_table[:, ["Symbol", feature_name]]
    dropmissing!(lookup_table)

    ElementFeatureDescriptor(
        feature_name,
        vector_length,
        nbins,
        logspaced,
        categorical,
        lookup_table
    )
end

# pretty printing, short version
Base.show(io::IO, af::ElementFeatureDescriptor) = print(io, "ElementFeature $(af.name)")

# pretty printing, long version
function Base.show(io::IO, ::MIME"text/plain", af::ElementFeatureDescriptor)
    st = "ElementFeature $(af.name):\n   categorical: $(af.categorical)\n   encoded length: $(af.length)"
    print(io, st)
end

# TODO: add way to get range/list of possible values for feature...

encodable_elements(f::ElementFeatureDescriptor) = f.lookup_table[:, :Symbol]

function encodable_elements(feature_name::String, lookup_table::DataFrame = atom_data_df)
    info = lookup_table[:, [Symbol(feature_name), :Symbol]]
    return info[
        findall(x -> !ismissing(x), getproperty(info, Symbol(feature_name))),
        :Symbol,
    ]
end

function (f::ElementFeatureDescriptor)(a::AbstractAtoms)
    @assert all([el in encodable_elements(f) for el in a.elements]) "Feature $(f.name) cannot encode some element(s) in this structure!"
    reduce(
        hcat,
            map(
                e -> onehot_lookup_encoder(
                    e,
                    f.name,
                    f.lookup_table;
                    nbins = f.nbins,
                    logspaced = f.logspaced,
                    categorical = f.categorical,
            ),
            a.elements,
        )
    )
end

#=
    decode_f =
        encoded_feature -> onecold_decoder(
            encoded_feature,
            feature_name,
            lookup_table;
            nbins = nbins,
            logspaced = logspaced,
            categorical = categorical,
        )
=#

function decode(f::ElementFeatureDescriptor, encoded_feature)
    onecold_decoder(encoded_feature, f.name, f.lookup_table;
                    nbins = f.nbins, logspaced = f.logspaced, categorical = f.categorical)
end


