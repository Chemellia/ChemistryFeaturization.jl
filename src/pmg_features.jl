#using PyCall
using JSON
using CSV
using DataFrames
include("atomfeat.jl")

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

"Little helper function to check that the logspace vector/boolean is appropriate and convert it to a vector as needed."
function get_logspaced_vec(vec::Vector, num_features::Integer)
    if vec==false # default behavior
        logspaced_vec = [false for i in 1:num_features]
    elseif vec==true
        logspaced_vec = [true for i in 1:num_features]
    elseif length(vec) == num_features # specified properly
        logspaced_vec = vec
    elseif length(vec) < num_features
        println("logspaced vector too short. Padding end with falses.")
        logspaced_vec = hcat(vec, [false for i in 1:num_features-size(vec,1)])
    elseif size(vec, 1) > num_features
        println("logspaced vector too long. Cutting off at appropriate length.")
        logspaced_vec = vec[1:num_features]
    end
    return logspaced_vec
end

"""
Function to build a featurization given vectors of metadata.
"""
function build_atom_feats(feature_names::Vector{Symbol}, nbins::Vector{<:Integer}=default_nbins*ones(Int64, size(feature_names,1)), logspaced=false)
    num_features = length(feature_names)

    # figure out spacing for each feature
    logspaced_vec = get_logspaced_vec(logspaced, num_features)

    # build AtomFeat objects
    # TODO: make this a preallocation since we know the length
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
make_feature_vectors(feature_names::Vector{Symbol}, nbins::Vector{<:Integer}=default_nbins*ones(Int64, size(feature_names,1)), logspaced=false) = make_feature_vectors(build_atom_feats(feature_names, nbins, logspaced))

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
    nbins = [length(f) for f in features]
    # First, check that the featurization is valid
    if !(vec_valid(vec, nbins))
        error("Vector is invalid!")
    else
        chunks = chunk_vec(vec, nbins)

        # make dict from feature names to corresponding chunks in vector
        fea_chunks = Dict(zip(features, chunks))

        return Dict(feature=>onecold_bins(f, fea_chunks[f]) for f in features)
    end
end

# alternate call signature
decode_feature_vector(vec::Vector{<:Real}, feature_names::Vector{String}, nbins::Vector{<:Integer}, logspaced=false) = decode_feature_vector(vec, build_atom_feats(feature_names, nbins, logspaced))
