#=
This module houses the built-in feature values for a variety of non-contextual atom features and also some convenience functions for constructing them easily.
=#

using DataFrames
using ..ChemistryFeaturization.Data: atom_data_df, feature_info

# export things
export avail_feature_names, categorical_feature_names, categorical_feature_vals, continuous_feature_names
import ..ChemistryFeaturization: default_log, default_categorical, get_bins

# read in features...
const categorical_feature_names = feature_info["categorical"]

const continuous_feature_names = feature_info["continuous"]
const avail_feature_names =
    cat(categorical_feature_names, continuous_feature_names; dims = 1)

"Compute the minimum and maximum possible values of a feature."
function fea_minmax(feature_name::String, lookup_table::DataFrame = atom_data_df)
    @assert feature_name in names(lookup_table) "Feature $feature_name isn't in the lookup table!"
    return [
        f(skipmissing(lookup_table[:, Symbol(feature_name)])) for f in [minimum, maximum]
    ]
end

# convenience dispatches of OHOC util functions to work directly with lookup table...
default_log(
    feature_name::String,
    lookup_table::DataFrame = atom_data_df;
    threshold_oom::Real = 2,
) = default_log(fea_minmax(feature_name, lookup_table)...; threshold_oom=threshold_oom)

function default_categorical(feature_name::String, lookup_table::DataFrame = atom_data_df; threshold_length = 5)
    local categorical
    if feature_name in avail_feature_names
        if feature_name in categorical_feature_names
            categorical = true
        else
            categorical = false
        end
    else
        @assert feature_name in colnames && "Symbol" in colnames "Your lookup table must have a column called :Symbol and one with the same name as your feature ($(feature_name)) to be usable!"
        categorical = default_categorical(lookup_table[:,Symbol(feature_name)], threshold_length=threshold_length)
    end
    return categorical
end

function get_bins(
    feature_name::String,
    lookup_table::DataFrame = atom_data_df;
    nbins::Integer = 10,
    threshold_oom=2,
    threshold_length=5,
    logspaced::Bool = default_log(feature_name, lookup_table; threshold_oom=threshold_oom),
    categorical::Bool = default_categorical(feature_name, lookup_table; threshold_length=threshold_length),
)
    colnames = names(lookup_table)
    @assert feature_name in colnames && "Symbol" in colnames "Your lookup table must have a column called :Symbol and one with the same name as your feature ($(feature_name)) to be usable!"

    possible_vals = unique(skipmissing(lookup_table[:,Symbol(feature_name)]))

    get_bins(possible_vals, nbins=nbins, logspaced=logspaced, categorical=categorical)
end

# this may just not be needed...
# "Little helper function to check that the logspace/categorical vector/boolean is appropriate and convert it to a vector as needed."
# function get_param_vec(vec, num_features::Integer; pad_val = false)
#     if !(typeof(vec) <: Vector)
#         output_vec = [vec for i = 1:num_features]
#     elseif length(vec) == num_features # specified properly
#         output_vec = vec
#     elseif length(vec) < num_features
#         @info "Parameter vector too short. Padding end with $pad_val."
#         output_vec = vcat(vec, [pad_val for i = 1:num_features-size(vec, 1)])
#     elseif size(vec, 1) > num_features
#         @info "Parameter vector too long. Cutting off at appropriate length."
#         output_vec = vec[1:num_features]
#     end
#     return output_vec
# end
