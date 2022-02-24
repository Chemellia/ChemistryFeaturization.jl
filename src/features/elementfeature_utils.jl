#=
This module houses the built-in feature values for a variety of non-contextual atom features and also some convenience functions for constructing them easily.
=#

using DataFrames
using ..ChemistryFeaturization.Data: element_data_df, elementfeature_info

# export things
export elementfeature_names,
    categorical_elementfeature_names, continuous_elementfeature_names
import ..ChemistryFeaturization: default_log, default_categorical, get_bins

# read in features...
const categorical_elementfeature_names = elementfeature_info["categorical"]

const continuous_elementfeature_names = elementfeature_info["continuous"]
const elementfeature_names =
    cat(categorical_elementfeature_names, continuous_elementfeature_names; dims = 1)

"""
    fea_minmax(feature_name, lookup_table)

Compute the minimum and maximum possible values of an `ElementFeatureDescriptor`, given a(n optional) lookup table.
"""
function fea_minmax(feature_name::String, lookup_table::DataFrame = element_data_df)
    @assert feature_name in names(lookup_table) "Feature $feature_name isn't in the lookup table!"
    return [
        f(skipmissing(lookup_table[:, Symbol(feature_name)])) for f in [minimum, maximum]
    ]
end

"""
    default_log(feature_name, lookup_table=element_data_df; threshold_oom=2)

Determine whether an element feature should be encoded by a `OneHotOneCold` codec with logarithmically spaced bins.
"""
default_log(
    feature_name::String,
    lookup_table::DataFrame = element_data_df;
    threshold_oom::Real = 2,
) = default_log(fea_minmax(feature_name, lookup_table)...; threshold_oom = threshold_oom)

"""
    default_categorical(feature_name, lookup_table=element_data_df; threshold_length=5)

Determine whether an element feature should be treated as categorical- or continuous-valued.
"""
function default_categorical(
    feature_name::String,
    lookup_table::DataFrame = element_data_df;
    threshold_length = 5,
)
    local categorical
    if feature_name in elementfeature_names
        if feature_name in categorical_elementfeature_names
            categorical = true
        else
            categorical = false
        end
    else
        colnames = names(lookup_table)
        @assert feature_name in colnames && "Symbol" in colnames "Your lookup table must have a column called :Symbol and one with the same name as your feature ($(feature_name)) to be usable!"
        categorical =
            default_categorical(lookup_table[:, Symbol(feature_name)], threshold_length)
    end
    return categorical
end

"""
    get_bins(feature_name, lookup_table=element_data_df; nbins=10, threshold_oom=2, threshold_length=5, logspaced, categorical)

Compute list of bins for an element feature.
"""
function get_bins(
    feature_name::String,
    lookup_table::DataFrame = element_data_df;
    nbins::Integer = 10,
    threshold_oom = 2,
    threshold_length = 5,
    logspaced::Bool = default_log(
        feature_name,
        lookup_table;
        threshold_oom = threshold_oom,
    ),
    categorical::Bool = default_categorical(
        feature_name,
        lookup_table;
        threshold_length = threshold_length,
    ),
)
    colnames = names(lookup_table)
    @assert feature_name in colnames && "Symbol" in colnames "Your lookup table must have a column called :Symbol and one with the same name as your feature ($(feature_name)) to be usable!"

    possible_vals = unique(skipmissing(lookup_table[:, Symbol(feature_name)]))

    get_bins(possible_vals, nbins = nbins, logspaced = logspaced, categorical = categorical)
end

"Little helper function to check that the logspace/categorical vector/boolean is appropriate and convert it to a vector as needed."
function get_param_vec(vec, num_features::Integer; pad_val = false)
    if !(typeof(vec) <: Vector)
        output_vec = [vec for i = 1:num_features]
    elseif length(vec) == num_features # specified properly
        output_vec = vec
    elseif length(vec) < num_features
        @info "Parameter vector too short. Padding end with $pad_val."
        output_vec = vcat(vec, [pad_val for i = 1:num_features-size(vec, 1)])
    elseif size(vec, 1) > num_features
        @info "Parameter vector too long. Cutting off at appropriate length."
        output_vec = vec[1:num_features]
    end
    return output_vec
end
