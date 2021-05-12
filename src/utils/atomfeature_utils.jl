#=
This module houses the built-in feature values for a variety of non-contextual atom features and also some convenience functions for constructing them easily.
=#
module AtomFeatureUtils

using DataFrames
using CSV
using JSON
using Flux: onecold

# export things
export default_nbins, oom_threshold_log
export atom_data_df, avail_feature_names
export categorical_feature_names, categorical_feature_vals, continuous_feature_names
export default_log, fea_minmax, default_categorical
export get_bins, build_onehot_vec
export onehot_lookup_encoder, onecold_decoder, encodable_elements

# default number of bins for continuous features, if unspecified
const default_nbins = 10
# if values of a feature span more than this many orders of magnitude, log-space it by default (open to better names for this...)
const oom_threshold_log = 2

# read in features...
atom_data_path = joinpath(@__DIR__, "data", "pymatgen_atom_data.csv")
const atom_data_df = DataFrame(CSV.File(atom_data_path))
feature_info_path = joinpath(@__DIR__, "data", "feature_info.json")
const feature_info = JSON.parsefile(feature_info_path)

const categorical_feature_names = feature_info["categorical"]
const categorical_feature_vals = Dict(
    fea => sort(collect(Set(skipmissing(atom_data_df[:, fea])))) for
    fea in categorical_feature_names
)
# but I want blocks to be in my order
categorical_feature_vals["Block"] = ["s", "p", "d", "f"]
const continuous_feature_names = feature_info["continuous"]
const avail_feature_names =
    cat(categorical_feature_names, continuous_feature_names; dims = 1)

# helper function
function fea_minmax(feature_name::String, lookup_table::DataFrame = atom_data_df)
    @assert feature_name in names(lookup_table) "Feature $feature_name isn't in the lookup table!"
    return [
        f(skipmissing(lookup_table[:, Symbol(feature_name)])) for f in [minimum, maximum]
    ]
end

# helper function
function default_log(
    feature_name::String,
    lookup_table::DataFrame = atom_data_df;
    threshold::Real = oom_threshold_log,
)
    min_val, max_val = fea_minmax(feature_name, lookup_table)
    local log
    if typeof(min_val) <: Number
        signs = sign.(fea_minmax(feature_name, lookup_table))
        same_sign = all(x -> x == signs[1], signs)
        if same_sign
            oom_arg = sign(min_val) < 0 ? min_val / max_val : max_val / min_val
            oom = log10(oom_arg)
            log = oom > oom_threshold_log
        else
            log = false
        end
    else
        log = false
    end
    return log
end

# helper function - if no info, be categorical for non-numbers and noncategorical for numbers
function default_categorical(feature_name::String, lookup_table::DataFrame = atom_data_df)
    local categorical
    if feature_name in avail_feature_names
        if feature_name in categorical_feature_names
            categorical = true
        else
            categorical = false
        end
    else
        feature_type = eltype(skipmissing(lookup_table[:, Symbol(feature_name)]))
        if feature_type <: Number
            categorical = false
        else
            categorical = true
        end
    end
    return categorical
end

# helper function for encoder and decoder...
function get_bins(
    feature_name::String;
    lookup_table::DataFrame = atom_data_df,
    nbins::Integer = default_nbins,
    logspaced::Bool = default_log(feature_name, lookup_table),
    categorical::Bool = default_categorical(feature_name, lookup_table),
)
    local bins, min_val, max_val

    if categorical
        bins = unique(lookup_table[:, Symbol(feature_name)])
    else
        min_val, max_val = fea_minmax(feature_name, lookup_table)

        if isapprox(min_val, max_val)
            @warn "It looks like the minimum and maximum possible values of $feature_name are approximately equal. This could cause numerical issues with binning, and also this feature is likely uninformative. Perhaps reconsider if it needs to be included?"
        end

        if logspaced
            @assert all(x -> x == x[1], [min_val, max_val]) "I don't know how to do a logarithmically spaced feature whose value can be zero! :("
            if sign(min_val) > 0
                bins = 10 .^ range(log10(min_val), log10(max_val), length = nbins + 1)
            else
                bins =
                    -1 .* (
                        10 .^
                        range(log10(abs(min_val)), log10(abs(max_val)), length = nbins + 1)
                    )
            end
        else
            bins = range(min_val, max_val, length = nbins + 1)
        end
    end
    return bins
end


# another helper function
function build_onehot_vec(val, bins, categorical)
    local bin_index, onehot_vec
    if categorical
        onehot_vec = [0.0 for i = 1:length(bins)]
        bin_index = findfirst(isequal(val), bins)
    else
        onehot_vec = [0.0 for i = 1:(length(bins)-1)]
        bin_index = searchsorted(bins, val).stop
        if bin_index == length(bins) # got the max value
            bin_index = bin_index - 1
        elseif isapprox(val, bins[1]) # sometimes will get 0 if this doesn't get checked
            bin_index = 1
        end
    end
    onehot_vec[bin_index] = 1.0
    return onehot_vec
end

# docstring
function onehot_lookup_encoder(
    el::String,
    feature_name::String;
    lookup_table::DataFrame = atom_data_df,
    nbins::Integer = default_nbins,
    logspaced::Bool = default_log(feature_name, lookup_table),
    categorical::Bool = default_categorical(feature_name, lookup_table),
)
    colnames = names(lookup_table)
    @assert (feature_name in colnames) && ("Symbol" in colnames) "Your lookup table must have a column called :Symbol and one with the same name as your feature to be usable!"
    @assert (feature_name in colnames) && ("Symbol" in colnames) "Your lookup table must have a column called :Symbol and one with the same name as your feature to be usable!"

    feature_vals = lookup_table[:, [:Symbol, Symbol(feature_name)]]

    @assert el in feature_vals.Symbol "Element $el is not in the database! :("

    bins = get_bins(
        feature_name;
        nbins = nbins,
        logspaced = logspaced,
        categorical = categorical,
        lookup_table = lookup_table,
    )

    # pull value of feature for this element
    val = getproperty(feature_vals[feature_vals.Symbol.==el, :][1, :], Symbol(feature_name))
    build_onehot_vec(val, bins, categorical)
end

# docstring
function onecold_decoder(
    encoded,
    feature_name::String;
    nbins::Integer = default_nbins,
    lookup_table::DataFrame = atom_data_df,
    logspaced::Bool = default_log(feature_name, lookup_table),
    categorical::Bool = default_categorical(feature_name, lookup_table),
)
    colnames = names(lookup_table)
    @assert feature_name in colnames && "Symbol" in colnames "Your lookup table must have a column called :Symbol and one with the same name as your feature to be usable!"

    bins = get_bins(
        feature_name;
        nbins = nbins,
        logspaced = logspaced,
        categorical = categorical,
        lookup_table = lookup_table,
    )
    local decoded
    if categorical # return value
        decoded = onecold(encoded, bins)
    else # return bounds (TODO: should this be a tuple or a vector..? I like tuple for distinguishing from encoded vectors, but it doesn't play so nice with broadcasting...)
        decoded = (onecold(encoded, bins[1:end-1]), onecold(encoded, bins[2:end]))
    end
    return decoded
end

# docstring
function encodable_elements(feature_name::String, lookup_table::DataFrame = atom_data_df)
    info = lookup_table[:, [Symbol(feature_name), :Symbol]]
    return info[
        findall(x -> !ismissing(x), getproperty(info, Symbol(feature_name))),
        :Symbol,
    ]
end

end