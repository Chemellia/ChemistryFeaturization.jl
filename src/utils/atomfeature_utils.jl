#=
This module exports the built-in feature values for a variety of non-contextual atom features and also some convenience functions for constructing them easily.
=#
module AtomFeatureUtils

using DataFrames
using CSV
using JSON
using ...ChemistryFeaturization: AtomFeature

# export new constructors
export AtomFeature

# default number of bins for continuous features, if unspecified
const default_nbins = 10

# read in features...
atom_data_path = joinpath(@__DIR__, "data", "pymatgen_atom_data.csv")
const atom_data_df = DataFrame(CSV.File(atom_data_path))
feature_info_path = joinpath(@__DIR__, "data", "feature_info.json")
const feature_info = JSON.parsefile(feature_info_path)

const categorical_feature_names = feature_info["categorical"]
const categorical_feature_vals = Dict(fea=>sort(collect(Set(skipmissing(atom_data_df[:, fea])))) for fea in categorical_feature_names)
# but I want blocks to be in my order
categorical_feature_vals["Block"] = ["s", "p", "d", "f"]
const continuous_feature_names = feature_info["continuous"]
const not_features = feature_info["not_features"] # atomic name, symbol
const avail_feature_names = cat(categorical_feature_names, continuous_feature_names; dims=1)

# compile min and max values of each feature...
const fea_minmax = Dict{String, Tuple{Real, Real}}()
for feature in avail_feature_names
    if !(feature in categorical_feature_names)
        minval = minimum(skipmissing(atom_data_df[:, feature]))
        maxval = maximum(skipmissing(atom_data_df[:, feature]))
        fea_minmax[feature] = (minval, maxval)
    end
end

#=
 goal: AtomFeat constructor that can just take in the name, or optionally other arguments for numbers of bins, etc. and construct encode_f and decode_f

 as in...
X = AtomFeat("X")
X = AtomFeat("X", numbins=12)
X = AtomFeat("X", logspaced=true)
=#

# docstring
function onehot_lookup_encoder(feature_name; nbins=default_nbins, logspaced=false)
    # TODO: write this
    # it can determine if it's categorical or not
    # should return a function
end

# docstring
function onecold_lookup_decoder(feature_name; nbins=default_nbins, logspaced=false)
    # TODO: write this
    # like above, returns a function
end

# TODO: add optional user-provided lookup table (will need to extend continuous/categorical feature name/val lists, make local versions and reference those instead)
function AtomFeature(feature_name; nbins=default_nbins, logspaced=false)
    @assert feature_name in continuous_feature_names || feature_name in categorical_feature_names "Cannot automatically build AtomFeat for $feature_name; I can't find it in a lookup table!"

    categorical = feature_name in categorical_feature_names
    if categorical
        length = length(categorical_feature_vals[feature_name])
    else
        length = nbins
    end
    encode_f = onehot_lookup_encoder(feature_name; nbins=nbins, logspaced=logspaced)
    decode_f = onecold_lookup_decoder(feature_name; nbins=nbins, logspaced=logspaced)
    # TODO: get type parameters, or decide we don't need them
    AtomFeat(feature_name, encode_f, decode_f, categorical, false, length)
end

end