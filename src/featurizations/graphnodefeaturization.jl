export GraphNodeFeaturization
export encodable_elements, decode, chunk_vec

"""
    GraphNodeFeaturization(atom_features::Vector{AtomFeature})
    GraphNodeFeaturization(feature_names, lookup_table, nbins, logspaced, categorical)

A featurization for AtomGraph objects that encodes features associated with each node. Contains a collection of `AtomFeature` objects, and can be initialized by passing those, or by passing parameters for constructing them.

## Notes
The "convenience constructor" that builds the AtomFeature objects for you only supports noncontextual AtomFeatures (i.e. those that can be encoded solely based on elemental identity and hence tabulated in a lookup table). If you have contextual features that require non-lookup-table encoding functions, you currently must build those features yourself and use the default constructor.

## Required Arguments
- `feature_names::Vector{String}`: Names of each feature
- `lookup_table::DataFrame` (optional): source of data for features, if not included in built-in `atom_data_df`

## Keyword Arguments
- `nbins::Union{Vector{Integer},Integer}`: Number of bins to use for one-hot encoding of continuous-valued features. Will be ignored for categorical features.
- `logspaced::Union{Vector{Bool},Bool}`: Whether to logarithmically space the bins
- `categorical::Union{Vector{Bool},Bool}`: Whether each feature is categorical or continous-valued.
"""
struct GraphNodeFeaturization <: AbstractFeaturization
    atom_features::Vector{AtomFeature}
end

function GraphNodeFeaturization(
    feature_names::Vector{String};
    lookup_table::Union{DataFrame,Nothing} = nothing,
    nbins::Union{Vector{Integer},Integer,Nothing} = nothing,
    logspaced::Union{Vector{Bool},Bool,Nothing} = nothing,
    categorical::Union{Vector{Bool},Bool,Nothing} = nothing,
)
    num_features = length(feature_names)
    local lookup_table_here , logspaced_here , categorical_here , nbins_here
    if isnothing(lookup_table)
        lookup_table_here = atom_data_df
    else
        # need to merge them in case some data is in one place and some the other
        lookup_table_here = outerjoin(atom_data_df, lookup_table, on = :Symbol)
    end

    if isnothing(logspaced)
        logspaced_here = map(fn -> default_log(fn, lookup_table_here), feature_names)
    else
        logspaced_here = get_param_vec(logspaced, num_features)
    end

    if isnothing(categorical)
        categorical_here =
            map(fn -> default_categorical(fn, lookup_table_here), feature_names)
    else
        categorical_here = get_param_vec(categorical, num_features)
    end

    if isnothing(nbins)
        nbins_here = [default_nbins for i = 1:num_features]
    else
        nbins_here = get_param_vec(nbins, num_features, pad_val = default_nbins)
    end

    afs = map(zip(feature_names, nbins_here, logspaced_here, categorical_here)) do args
        AtomFeature(
            args[1],
            lookup_table_here,
            nbins = args[2],
            logspaced = args[3],
            categorical = args[4],
        )
    end
    GraphNodeFeaturization(afs)
end

# TODO: function to compute total vector length?

encodable_elements(fzn::GraphNodeFeaturization) =
    intersect([f.encodable_elements for f in fzn.atom_features]...)

"""
    chunk_vec(vec, nbins)

Helper function that divides up an already-constructed feature vector into "chunks" (one for each feature) of lengths specified by the vector nbins.

Sum of nbins should be equal to the length of vec.

# Examples
```jldoctest
julia> chunk_vec([1,0,0,1,0], [3,2])
2-element Array{Array{Bool,1},1}:
 [1, 0, 0]
 [1, 0]
 ```
"""
function chunk_vec(vec::Vector{<:Real}, nbins::Vector{<:Integer})
    chunks = fill(Bool[], size(nbins, 1))
    @assert length(vec) == sum(nbins) "Total number of bins doesn't match length of feature vector."
    for i = 1:size(nbins, 1)
        if i == 1
            start_ind = 1
        else
            start_ind = sum(nbins[1:i-1]) + 1
        end
        chunks[i] = vec[start_ind:start_ind+nbins[i]-1]
    end
    return chunks
end

function decode(fzn::GraphNodeFeaturization, encoded::Matrix{<:Real})
    num_atoms = size(encoded, 2)
    nbins = [f.length for f in fzn.atom_features]
    local decoded = Dict{Integer,Dict{String,Any}}()
    for i = 1:num_atoms
        #println("atom $i")
        chunks = chunk_vec(encoded[:, i], nbins)
        decoded[i] = Dict{String,Any}()
        for (chunk, f) in zip(chunks, fzn.atom_features)
            #println("    $(f.name): $(decode(f, chunk))")
            decoded[i][f.name] = decode(f, chunk)
        end
    end
    return decoded
end
