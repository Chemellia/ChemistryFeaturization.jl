export default_log, default_categorical, get_bins

"""
    default_log(min_val, max_val; threshold = 2)
    default_log(possible_vals; threshold = 2)

Determine whether a `OneHotOneCold` codec used with a particular `FeatureDescriptor` should have logarithmically spaced bins. 

Operates by comparing the ratio of the maximum to minimum values to a specified order-of-magnitude threshold.
"""
function default_log(min_val::T, max_val::T; threshold_oom = 2) where {T}
    local log
    if T <: Number
        signs = sign.([min_val, max_val])
        same_sign = all(x -> x == signs[1], signs)
        if same_sign
            oom_arg = sign(min_val) < 0 ? min_val / max_val : max_val / min_val
            oom = log10(oom_arg)
            log = oom > threshold_oom
        else
            log = false
        end
    else
        log = false
    end
    return log
end

# adding this dispatch for uniformity with next one
function default_log(possible_vals::Vector{T}, threshold_oom = 2) where {T}
    if T <: Number
        default_log(
            minimum(possible_vals),
            maximum(possible_vals),
            threshold_oom = threshold_oom,
        )
    else
        return false
    end
end

"""
    default_categorical(possible_vals; threshold_length=5)

Determine if a feature should be treated as categorical or continuous-valued.

If the value type is not a number, always returns true. If it is numerical, returns true if the number of possible values is less than `threshold_length` and false otherwise.
"""
function default_categorical(possible_vals::Vector{T}, threshold_length = 5) where {T}
    num_vals = length(unique(possible_vals))
    if num_vals < threshold_length
        return true
    else
        if T <: Number
            return false
        else
            return true
        end
    end
end

"""
    get_bins(possible_vals; threshold_oom=2, threshold_length=5, nbins=10, logspaced, categorical)

Given a list of possible values, return a list of bins, making sensible default choices for the binning parameters.

See also: [`default_log`](@ref), [`default_categorical`](@ref)
"""
function get_bins(
    possible_vals::Vector;
    threshold_oom = 2,
    threshold_length = 5,
    nbins = 10,
    logspaced::Bool = default_log(possible_vals, threshold_oom),
    categorical::Bool = default_categorical(possible_vals, threshold_length),
)
    local bins
    if categorical
        bins = sort(unique(possible_vals))
    else
        min_val = minimum(possible_vals)
        max_val = maximum(possible_vals)

        if isapprox(min_val, max_val)
            @warn "It looks like the minimum and maximum possible values of $feature_name are approximately equal. This could cause numerical issues with binning, and also this feature is likely uninformative. Perhaps reconsider if it needs to be included?"
        end

        if logspaced
            @assert all(x -> sign(x) == sign(min_val), [min_val, max_val]) "I don't know how to do a logarithmically spaced feature whose value can be zero! :("
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
