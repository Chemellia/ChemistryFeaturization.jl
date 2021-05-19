using ..ChemistryFeaturization.Utils.AtomFeatureUtils
using DataFrames

"""
    AtomFeature(feature_name, encode_f, decode_f, categorical, contextual, length, encodable_elements)

Construct a feature object that encodes features associated with individual atoms. If `contextual==false`, then the encoding is done solely based on elemental identity, otherwise it depends on the neighborhood.

## Arguments
- `name::String`: the name of the feature
- `encode_f::Function`: a function that takes in <:AbstractAtoms and returns encoded values of this feature for the atoms in that structure
- `decode_f::Function`: inverse function to `encode_f`, takes in encoded feature and returns value (for categorical) or range of values (for continuous-valued) of the feature
- `categorical::Bool`: flag for whether the feature is categorical or continuous-valued
- `contextual::Bool`: flag for whether the feature's value depends on an atom's environment
- `length::Int`: length of encoded vector
- `encodable_elements::Vector{String}`: list of elemental symbols representing species that can be encoded by this feature
"""
struct AtomFeature <: AbstractFeature
    name::String
    encode_f::Function
    decode_f::Function
    categorical::Bool
    contextual::Bool
    length::Integer
    encodable_elements::Vector{String}
end

"""
    AtomFeature(feature_name, lookup_table = atom_data_df; nbins, logspaced, categorical)

Construct a non-contextual AtomFeature that encodes using a lookup table.

## Required Arguments
- `feature_name::String`: Name of the feature
- `lookup_table::DataFrame` (optional): if feature is not included in the built-in `atom_data_df`, provide a lookup table that includes its value for every element you want to encode it on

## Keyword Arguments
- `nbins::Integer`: Number of bins for one-hot encoding of continuous-valued features. Will be ignored for categorical features.
- `logspaced::Bool`: flag for whether bins should be logarithmically spaced
- `categorical::Bool`: flag for whether feature is categorical or continuous-valued

## Examples

"""
function AtomFeature(
    feature_name::String,
    lookup_table::DataFrame = atom_data_df;
    nbins::Integer = default_nbins,
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
        encodable_elements(
            feature_name,
            lookup_table,
        ),
    )
end

# pretty printing, short version
Base.show(io::IO, af::AtomFeature) = print(io, "AtomFeature $(af.name)")

# pretty printing, long version
function Base.show(io::IO, ::MIME"text/plain", af::AtomFeature)
    st = "AtomFeature $(af.name):\n   categorical: $(af.categorical)\n   contextual: $(af.contextual)\n   encoded length: $(af.length)"
    print(io, st)
end

"""
    feature_range(feature_name::String, lookup_table::DataFrame = atom_data_df)
    feature_range(feature::AtomFeature)

Return minimum and maximum values of the provided feature.
"""
function feature_range(feature_name::String, lookup_table::DataFrame = atom_data_df)
    @assert feature_name in names(lookup_table) "Feature $feature_name isn't in the lookup table!"
    return [
        f(skipmissing(lookup_table[:, Symbol(feature_name)])) for f in [minimum, maximum]
    ]
end

encodable_elements(f::AtomFeature) = f.encodable_elements

function encodable_elements(feature_name::String, lookup_table::DataFrame = atom_data_df)
    info = lookup_table[:, [Symbol(feature_name), :Symbol]]
    return info[
        findall(x -> !ismissing(x), getproperty(info, Symbol(feature_name))),
        :Symbol,
    ]
end

function (f::AtomFeature)(a::AbstractAtoms)
    @assert all([el in f.encodable_elements for el in a.elements]) "Feature $(f.name) cannot encode some element(s) in this structure!"
    f.encode_f(a)
end

# TODO: some Weave stuff needed here
