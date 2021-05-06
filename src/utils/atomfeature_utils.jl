#=
This module houses the built-in feature values for a variety of non-contextual atom features and also some convenience functions for constructing them easily.
=#
module AtomFeatureUtils

using DataFrames
using CSV
using JSON
using Flux: onecold

# export things
export default_nbins, atom_data_df, avail_feature_names, not_features
export categorical_feature_names, categorical_feature_vals 
export continuous_feature_names, fea_minmax
export get_bins, build_onehot_vec
export onehot_lookup_encoder, onecold_decoder

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

# helper function for encoder and decoder...
function get_bins(feature_name; nbins=default_nbins, logspaced=false)
    categorical = feature_name in categorical_feature_names
    local bins
    if categorical
        bins = categorical_feature_vals[feature_name]
    else
        min_val, max_val = fea_minmax[feature_name]
        if logspaced
            bins = 10 .^ range(log10(min_val), log10(max_val), length=nbins+1)
        else
            bins = range(min_val, max_val, length=nbins+1)
        end
    end
    return bins
end

# another helper function
function build_onehot_vec(val, bins, categorical)
    local bin_index, onehot_vec
    if categorical
        onehot_vec = [0.0 for i in 1:length(bins)]
        bin_index = findfirst(isequal(val), bins)
    else
        onehot_vec = [0.0 for i in 1:(length(bins)-1)]
        bin_index = searchsorted(bins, val).stop
        if bin_index == length(bins) # got the max value
            bin_index = bin_index - 1
        elseif isapprox(val, bins[1]) # sometimes will get 0 if this doesn't get checked
            bin_index=1
        end
    end
    onehot_vec[bin_index] = 1.0
    return onehot_vec
end

# docstring
function onehot_lookup_encoder(el::String, feature_name; nbins=default_nbins, logspaced=false)
    @assert feature_name in avail_feature_names "$feature_name is not a built-in feature, you'll have to write your own encoder function. Available built-in features are: $avail_feature_names"
    @assert el in atom_data_df.Symbol "Element $el is not in the database! :("

    feature_vals = atom_data_df[:, [:Symbol, Symbol(feature_name)]]
    categorical = feature_name in categorical_feature_names
    bins = get_bins(feature_name; nbins=nbins, logspaced=logspaced)

    # pull value of feature for this element
    val = getproperty(feature_vals[feature_vals.Symbol.==el,:][1,:], Symbol(feature_name))
    build_onehot_vec(val, bins, categorical)
end

# docstring
function onecold_decoder(encoded, feature_name; nbins=default_nbins, logspaced=false)
    @assert feature_name in avail_feature_names "$feature_name is not a built-in feature, you'll have to write your own decoder function. Available built-in features are: $avail_feature_names"

    bins = get_bins(feature_name; nbins=nbins, logspaced=logspaced)
    categorical = feature_name in categorical_feature_names

    if categorical # return value
        decoded = onecold(encoded, bins)
    else # return bounds
        decoded = (onecold(encoded, bins[1:end-1]), onecold(encoded, bins[2:end]))
    end
    return decoded
end

# TODO: add optional user-provided lookup table (will need to extend continuous/categorical feature name/val lists, make local versions and reference those instead)

end