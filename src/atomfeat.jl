using Flux: onecold
using JSON
using CSV
using DataFrames

# default number of bins for continuous features, if unspecified
const default_nbins = 10

# read in features...
atom_data_path = joinpath(dirname(pathof(ChemistryFeaturization)), "..", "data", "pymatgen_atom_data.csv")
const atom_data_df = DataFrame!(CSV.File(atom_data_path))
feature_info_path = joinpath(dirname(pathof(ChemistryFeaturization)), "..", "data", "feature_info.json")
const feature_info = JSON.parsefile(feature_info_path)

const categorical_feature_names = Symbol.(feature_info["categorical"])
const categorical_feature_vals = Dict(fea=>sort(collect(Set(skipmissing(atom_data_df[:, fea])))) for fea in categorical_feature_names)
# but I want blocks to be in my order
categorical_feature_vals[:Block] = ["s", "p", "d", "f"]
const continuous_feature_names = Symbol.(feature_info["continuous"])
const not_features = Symbol.(feature_info["not_features"]) # atomic name, symbol
const avail_feature_names = cat(categorical_feature_names, continuous_feature_names; dims=1)

# compile min and max values of each feature...
const fea_minmax = Dict{Symbol, Tuple{Real, Real}}()
for feature in avail_feature_names
    if !(feature in categorical_feature_names)
        minval = minimum(skipmissing(atom_data_df[:, feature]))
        maxval = maximum(skipmissing(atom_data_df[:, feature]))
        fea_minmax[feature] = (minval, maxval)
    end
end

# Type to store featurization metadata. An array of these specifies a featurization scheme for atoms.
# TODO for Sean, eventually: add PairFeat and BondFeat types (maybe in another file for tidiness)
# TODO, maybe: we could define AbstractFeat{T} that all of these would inherit from, but TBD if that would be useful or not
struct AtomFeat{T}
    name::Symbol # name of feature
    categorical::Bool # whether it's categorical (vs. numerical)
    num_bins::Integer # length of associated subvector
    logspaced::Bool # whether it's logspaced (most relevant for numerical)
    vals::Vector{T} # list of values (length equal to num_bins for categorical, num_bins+1 (specifying bin edges) for numerical)
    # basic standard constructor will check some things...
    function AtomFeat(name::Symbol, categorical::Bool, num_bins::Integer, logspaced::Bool, vals::Vector)
        T = eltype(vals)
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
            if !(eltype(vals)<:Real)
                error("Numerical features must have (real) numerical values...") # should figure out how to throw a TypeError propertly
            end
            val_list = sort(vals)
        end
        new{T}(name, categorical, num_bins, logspaced, val_list)
    end
end

# constructor for features where vals gets calculated rather than passed in, works for categorical too if their values are numbers
function AtomFeat(name::Symbol, categorical::Bool, num_bins::Integer, min_val::Real, max_val::Real, logspaced::Bool=false; T=Float32)
    categorical ? len = num_bins : len = num_bins + 1
    if logspaced
        vals = 10 .^ range(log10(min_val), log10(max_val), length=len)
    else
        vals = range(min_val, max_val, length=len)
    end
    AtomFeat(name, categorical, num_bins, logspaced, T.([v for v in vals]))
end

# constructor that will assume categorical features
AtomFeat(name::Symbol, vals::Vector) = AtomFeat(name, true, length(vals), false, vals)

# pretty printing, short form
Base.show(io::IO, f::AtomFeat{T}) where {T} = print(io, "$(f.name): AtomFeat{$T} with $(f.num_bins) bins")

# pretty printing, long form
Base.show(io::IO, ::MIME"text/plain", f::AtomFeat{T}) where{T} = print(io, "AtomFeat{$T}\n   name: $(f.name)\n   categorical: $(f.categorical)\n   length: $(f.num_bins)\n   logspaced: $(f.logspaced)\n   bins: $(f.vals)")

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
        elseif isapprox(val, f.vals[1])
            # sometimes it gives 0 otherwise
            bin_index = 1
        end
    end
    onehot_vec = [0.0 for i in 1:f.num_bins]
    onehot_vec[bin_index] = 1.0
    return onehot_vec
end

# TODO, maybe: add signatures for this and onecold version where you give the values rather than the feature object?

"""
    onecold_bins(feature, vec)

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
        decoded = decoded = (onecold(vec, f.vals[1:end-1]), onecold(vec, f.vals[2:end]))
    end
    return decoded
end

"Little helper function to check that the logspace vector/boolean is appropriate and convert it to a vector as needed."
function get_logspaced_vec(vec, num_features::Integer)
    if vec==false # default behavior
        logspaced_vec = [false for i in 1:num_features]
    elseif vec==true
        logspaced_vec = [true for i in 1:num_features]
    elseif length(vec) == num_features # specified properly
        logspaced_vec = vec
    elseif length(vec) < num_features
        println("logspaced vector too short. Padding end with falses.")
        logspaced_vec = vcat(vec, [false for i in 1:num_features-size(vec,1)])
    elseif size(vec, 1) > num_features
        println("logspaced vector too long. Cutting off at appropriate length.")
        logspaced_vec = vec[1:num_features]
    end
    return logspaced_vec
end

"""
Function to build a featurization given vectors of metadata.

Note that nbins will be ignored for categorical features.
"""
function build_atom_feats(feature_names::Vector{Symbol}; nbins::Vector{<:Integer}=default_nbins*ones(Int64, size(feature_names,1)), logspaced=false)
    num_features = length(feature_names)

    # figure out spacing for each feature
    logspaced_vec = get_logspaced_vec(logspaced, num_features)

    # build AtomFeat objects
    # TODO: could make this a preallocation since we know the length, probably not critical though
    feature_specs = AtomFeat[]
    #for i in 1:num_features
    for t in zip(feature_names, nbins, logspaced_vec)
        feature_name = t[1]
        nbins = t[2]
        logspaced = t[3]
        if feature_name in categorical_feature_names
            vals = categorical_feature_vals[Symbol(feature_name)]
            feature = AtomFeat(feature_name, vals)
        elseif feature_name in continuous_feature_names
            feature = AtomFeat(feature_name, false, nbins, fea_minmax[feature_name]..., logspaced)
        end
        push!(feature_specs, feature)
    end
    return feature_specs
end

"""
    make_feature_vectors(features)
    make_feature_vectors(feature_names)
    make_feature_vectors(feature_names, nbins)
    make_feature_vectors(feature_names, nbins, logspaced)


Make custom feature vectors, using specified features and numbers of bins. Can be called with an array of AtomFeat objects or by specifying names and other metadata and the objects will be built.

Note that in the latter case, bin numbers will be ignored for categorical features (block, group, and row), but features and nbins vectors should still be the same length. Optionally, feed in vector of booleans with trues at the index of any (continous valued) feature whose bins should be log spaced.

# Arguments
- `feature_names::Vector{String}`: list of features to be encoded
- `nbins::Vector{Integer}`: number of bins for each feature (in same order)
- `logspaced=false`: (single Bool or vector of them) whether or not to logarithmically space each feature

# Returns 
- a dictionary from element symbol => one-hot style feature vector, concatenated in order of feature list.
- A Vector of AtomFeat objects
"""
function make_feature_vectors(features::Vector{AtomFeat})
    feature_names = [f.name for f in features]
    usable_atom_data =  dropmissing(atom_data_df[:, cat(Symbol.(feature_names), not_features, dims=1)])
    num_dropped = size(atom_data_df)[1] - size(usable_atom_data)[1]
    if num_dropped != 0
        @info "$(num_dropped) elements were dropped so that all features are defined."
    end

    # dict from element symbol to feature vec of that element
    # (if we do any structure-specific features later, e.g. coordination or something,
    # this will have to iterate over every atom in the structure instead...)
    # (but possibly would just want to append those to the end anyway...)
    sym_featurevec = Dict{String, Vector{Float32}}()
    
    for i in 1:size(usable_atom_data,1)
        el = usable_atom_data.Symbol[i]
        featurevec = []
        # make onehot vector for each feature
        for f in features
            val = usable_atom_data[i, Symbol(f.name)]
            subvec = onehot_bins(f, val)
            append!(featurevec, subvec)
        end
        sym_featurevec[el] = featurevec
    end
    return sym_featurevec, features
end

# alternate call signature
make_feature_vectors(feature_names::Vector{Symbol}, nbins::Vector{<:Integer}=default_nbins*ones(Int64, size(feature_names,1)), logspaced=false) = make_feature_vectors(build_atom_feats(feature_names; nbins=nbins, logspaced=logspaced))

"""
    chunk_vec(vec, nbins)

Divide up an already-constructed feature vector into "chunks" (presumably one for each feature) of lengths specified by the vector nbins.

Sum of nbins should be equal to the length of vec.

# Examples
```jldoctest
julia> chunk_vec([1,0,0,1,0], [3,2])
2-element Array{Array{Bool,1},1}:
 [1, 0, 0]
 [1, 0]
 ```
"""
function chunk_vec(vec::Vector, nbins::Vector{<:Integer})
    chunks = fill(Bool[], size(nbins, 1))
    @assert length(vec)==sum(nbins) "Total number of bins doesn't match length of feature vector."
    for i in 1:size(nbins,1)
        if i==1
            start_ind = 1
        else
            start_ind = sum(nbins[1:i-1])+1
        end
        chunks[i] = vec[start_ind:start_ind+nbins[i]-1]
    end
    return chunks
end

"""
Check that each subvector in a featurization is a valid one-hot encoding (one true and otherwise all falses).
"""
function vec_valid(vec::Vector, nbins::Vector{<:Integer})
    result = true
    # are there the right number of total bins?
    try
        chunks = chunk_vec(vec, nbins)
    catch e
        result = false
    end
    if result # if we're still good up to now...
        chunks = chunk_vec(vec, nbins) # apparently chunks is local to the previous try block?
        # does each subvector have exactly one true value?
        for i in 1:length(nbins)
            subvec = chunks[i]
            if !(sum(subvec)==1)
                println("Subvector ", i, " is invalid.")
                result = false
            end
        end
    end
    return result
end

"""
Function to invert the binning process. Useful to check that it's working properly, or just to inspect properties once they've been encoded.

Need to feed in a feature vector as well as the lists of features and bin numbers that were used to encode it.
"""
function decode_feature_vector(vec::Vector, features::Vector{AtomFeat})
    nbins = [f.num_bins for f in features]
    # First, check that the featurization is valid
    if !(vec_valid(vec, nbins))
        error("Vector is invalid!")
    else
        chunks = chunk_vec(vec, nbins)

        # make dict from feature names to corresponding chunks in vector
        fea_chunks = Dict(zip(features, chunks))

        return Dict(f.name=>onecold_bins(f, fea_chunks[f]) for f in features)
    end
end

# alternate call signature
decode_feature_vector(vec::Vector{<:Real}, feature_names::Vector{String}, nbins::Vector{<:Integer}, logspaced=false) = decode_feature_vector(vec, build_atom_feats(feature_names, nbins, logspaced))