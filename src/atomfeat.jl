using Flux: onecold

# Type to store featurization metadata. An array of these specifies a featurization scheme for atoms.
# TODO for Sean, eventually: add PairFeat and BondFeat types
struct AtomFeat{T}
    name::Symbol
    categorical::Bool
    num_bins::Integer
    logspaced::Bool
    vals::Vector{T}
    # basic standard constructor will check some things...
    function AtomFeat(name::Symbol, categorical::Bool, num_bins::Integer, logspaced::Bool, vals::Vector)
        T = typeof(vals[1])
        if categorical
            if num_bins != length(vals)
                DimensionMismatch("Categorical features should have a number of values equal to the number of bins.")
            end
            if logspaced
                warn("A logarithmically spaced categorical feature doesn't make much sense...are you sure that's what you meant?")
            end
            val_list = vals
        else # numerical feature
            if num_bins + 1 != length(vals)
                DimensionMismatch("Numerical features should have values specifying bin edges. This vector is the wrong length!")
            end
            if !eltype(vals)<:Real)
                error("Numerical features must have (real) numerical values...") # should figure out how to throw a TypeError propertly
            end
            val_list = sort(vals)
        end
        new{T}(name, categorical, num_bins, logspaced, val_list)
    end
end

# constructor for features where vals gets calculated rather than passed in, works for categorical too if their values are numbers
function AtomFeat(name::Symbol, categorical::Bool, num_bins::Integer, min_val::Real, max_val::Real, logspaced::Bool=false)
    categorical ? len = num_bins : len = num_bins + 1
    if logspaced
        vals = 10 .^ range(log10(min_val), log10(max_val), length=len)
    else
        vals = range(min_val, max_val, length=len)
    end
    AtomFeat(name, categorical, num_bins, logspaced, [v for v in vals])
end

# constructor that will assume categorical features
AtomFeat(name::Symbol, vals::Vector) = AtomFeat(name, true, length(vals), false, vals)

"""
Create onehot style vector.

# Arguments
- `feature::AtomFeat`: feature being encoded
- `val`: value of feature to encode

See also: [onecold_bins](@ref), [onehot](@ref Flux.onehot)
"""
function onehot_bins(f::AtomFeat, val)
    if f.categorical
        bin_index = findfirst(isequal(val), f.vals)
    else
        bin_index = searchsorted(f.vals, val).stop
        if bin_index == length(f.vals) # got the max value
            bin_index = bin_index - 1
        end
    end
    onehot_vec = [0.0 for i in 1:f.num_bins]
    onehot_vec[bin_index] = 1.0
    return onehot_vec
end

"""
    onecold_bins(feature, vec, bins)

Inverse function to onehot_bins, decodes a vector corresponding to a given feature. For a categorical feature, returns a value, for a continuous one, a range.

# Arguments
- `feature::AtomFeat`: feature being encoded
- `vec::Vector`: vector (such as produced by [onehot_bins](@ref))

See also: [onehot_bins](@ref), [onecold](@ref Flux.onecold)

"""
function onecold_bins(f::AtomFeat, vec::Vector)
    if f.categorical
        decoded = onecold(vec, f.vals)
    else
        decoded = decoded = (onecold(vec, f.vals[1:end-1]), onecold(vec, vals[2:end]))
    end
    return decoded
end
