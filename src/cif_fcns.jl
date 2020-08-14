using PyCall
using JSON
using CSV
using DataFrames
using Flux: onecold

# default number of bins for numerical features, if unspecified
const default_nbins = 10

# read in features...
feature_info = JSON.parsefile("../data/feature_info.json")
global categorical_features = feature_info["categorical"]
categorical_feature_vals = Dict(fea=>sort(collect(Set(skipmissing(atom_data_df[:, fea])))) for fea in categorical_features)
# but I want blocks to be in my order
categorical_feature_vals["Block"] = ["s", "p", "d", "f"]
global numerical_features = feature_info["numerical"]
global not_features = feature_info["not_features"] # atomic name, symbol

avail_features = cat(categorical_features, numerical_features; dims=1)

atom_data_df = DataFrame!(CSV.File("../data/pymatgen_atom_data.csv"))

# compile min and max values of each feature...
global fea_minmax = Dict()
for feature in avail_features
    if !(feature in categorical_features)
        minval = minimum(skipmissing(atom_data_df[:, Symbol(feature)]))
        maxval = maximum(skipmissing(atom_data_df[:, Symbol(feature)]))
        fea_minmax[feature] = (minval, maxval)
    end
end

"Get bins for a given feature, intelligently handling categorical vs. continuous feature values. In the former case, returns the categories. In the later, returns bin edges."
function get_bins(feature; nbins=default_nbins, logspaced=false)
    if feature in categorical_features
        bins = categorical_feature_vals[feature]
    else
        min_val = fea_minmax[feature][1]
        max_val = fea_minmax[feature][2]
        if logspaced
            bins = 10 .^range(log10(min_val), log10(max_val), length=nbins+1)
        else
            bins = range(min_val, max_val, length=nbins+1)
        end
    end
    return bins
end

"Find which bin index a value sits in, intelligently handling both categorical and continous feature values."
function which_bin(feature, val, bins=get_bins(feature))
    if feature in categorical_features
        bin_index = findfirst(isequal(val), bins)
    else
        bin_index = searchsorted(bins, val).stop
        if bin_index == length(bins) # got the max value
            bin_index = bin_index-1
        end
    end
    return bin_index
end

"""
    onehot_bins(feature, val)
    onehot_bins(feature, val, bins)

Create onehot style vector, handling both categorical and continuous features.

# Arguments
- `feature::String`: feature being encoded
- `val`: (Float or String) value of feature to encode
- `bins::Array`: categorical values or bin edges if continouus feature

# Examples
```jldoctest
julia> onehot_bins("cont_feature", 3, [0,2,4,6])
3-element Array{Bool,1}:
 0
 1
 0

julia> onehot_bins("cat_feature", 3, [1,2,3])
3-element Array{Bool,1}:
 0
 0
 1
```

See also: [onecold_bins](@ref), [onehot](@ref Flux.onehot)

"""
function onehot_bins(feature, val, bins=get_bins(feature))
    if feature in categorical_features
        len = length(bins)
    else
        len = length(bins)-1
    end
    onehot_vec = [0.0 for i in 1:len]
    onehot_vec[which_bin(feature, val, bins)] = 1.0
    return onehot_vec
end

"""
    onecold_bins(feature, vec, bins)

Inverse function to onehot_bins, decodes a vector corresponding to a given feature, given the bins that were used to encode it.

# Arguments
- `feature::String`: feature being encoded
- `vec::Array`: vector (such as produced by [onehot_bins](@ref))
- `bins::Array`: categorical values or bin edges if continouus feature

# Examples
```jldoctest
julia> onecold_bins("cont_feature", [0,1,0], [0,2,4,6])
(2,4)

julia> onecold_bins("cat_feature", [0,0,1], [1,2,3])
3
```

See also: [onehot_bins](@ref), [onecold](@ref Flux.onecold)

"""
function onecold_bins(feature, vec, bins)
    if feature in categorical_features
        # return value
        decoded = onecold(vec, bins)
    else
        # return range of values
        decoded = (onecold(vec, bins[1:end-1]), onecold(vec, bins[2:end]))
    end
    return decoded
end


"Little helper function to check that the logspace vector/boolean is appropriate and convert it to a vector as needed."
function get_logspaced_vec(vec, num_features)
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
    make_feature_vectors(features)
    make_feature_vectors(features, nbins)
    make_feature_vectors(features, nbins, logspaced)

Make custom feature vectors, using specified features and numbers of bins. Note that bin numbers will be ignored for categorical features (block, group, and row), but features and nbins vectors should still be the same length (there's probably a more elegant way to handle that).

Optionally, feed in vector of booleans with trues at the index of any (continous valued) feature whose bins should be log spaced.

# Arguments
- `features::Array{String,1}`: list of features to be encoded
- `nbins::Array{Integer,1}`: number of bins for each feature (in same order)
- `logspaced::Array{Bool,1}=false`: whether or not to logarithmically space each feature

Returns a dictionary from element symbol => one-hot style feature vector, concatenated in order of feature list.
"""
function make_feature_vectors(features, nbins=default_nbins*ones(Int64, size(features,1)), logspaced=false)
    num_features = size(features,1)

    # figure out spacing for each feature
    logspaced_vec = get_logspaced_vec(logspaced, num_features)

    # make dict from feature name to bins
    features_bins = Dict(features[i] => get_bins(features[i]; nbins=Int64(nbins[i]), logspaced=logspaced_vec[i]) for i in 1:num_features)

    # dict from feature name to number of bins for that feature
    features_nbins = Dict(zip(features, nbins))

    # dict from element symbol to feature vec of that element
    # (if we do any structure-specific features later, e.g. coordination or something,
    # this will have to iterate over every atom in the structure instead...)
    # (but possibly would just want to append those to the end anyway...)
    sym_featurevec = Dict{String, Array{Float32,1}}()
    usable_atom_data =  dropmissing(atom_data_df[:, cat(features, not_features, dims=1)])
    orig_len = size(atom_data_df)[1]
    new_len = size(usable_atom_data)[1]
    diff = orig_len - new_len
    if orig_len != new_len # had to drop some atoms
        @info "$(diff) elements were dropped so that all features are defined."
    end
    for i in 1:size(usable_atom_data,1)
        el = usable_atom_data.Symbol[i]
        featurevec = []
        # make onehot vector for each feature
        for feature in features
            feature_val = usable_atom_data[i, Symbol(feature)]
            subvec = onehot_bins(feature, feature_val, get_bins(feature; nbins=features_nbins[feature]))
            append!(featurevec, subvec)
        end
        sym_featurevec[el] = featurevec # need transpose because of how graphcon works
    end
    return sym_featurevec
end

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
function chunk_vec(vec, nbins)
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
function vec_valid(vec, nbins)
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
function decode_feature_vector(vec, features, nbins, logspaced=false)
    # First, check that the featurization is valid
    if !(vec_valid(vec, nbins))
        println("Vector is invalid!")
    else
        chunks = chunk_vec(vec, nbins)
        # make dict from features to corresponding chunks
        fea_chunks = Dict(zip(features, chunks))

        # and one from feature to bin bounds for this vector
        num_features = length(features)
        logspaced_vec = get_logspaced_vec(logspaced, num_features)
        fea_bins = Dict(features[i]=>get_bins(features[i]; nbins=nbins[i], logspaced=logspaced_vec[i]) for i in 1:num_features)

        return Dict(feature=>onecold_bins(feature, fea_chunks[feature], fea_bins[feature]) for feature in features)
    end
end


## now come the graph-building functions

# options for decay of bond weights with distance...
inverse_square(x) = x^-2.0
exp_decay(x) = exp(-x)

"""
Function to build graph from a CIF file of a crystal structure. Returns a tuple of the adjacency matrix and the list of elements in order of the graph nodes.

# Arguments
- `cif_path::String`: path to CIF file
- `use_voronoi::bool`: if true, use Voronoi method for neighbor lists, if false use cutoff method

    (The rest of these parameters are only used if use_voronoi is false)

- `radius::Float=8.0`: cutoff radius for atoms to be considered neighbors (in angstroms)
- `max_num_nbr::Integer=12`: maximum number of neighbors to include (even if more fall within cutoff radius)
- `dist_decay_func`: function (e.g. inverse_square or exp_decay) to determine falloff of graph edge weights with neighbor distance
"""
function build_graph(cif_path; use_voronoi=true, radius=8.0, max_num_nbr=12, dist_decay_func=inverse_square, normalize=true)
    s = pyimport("pymatgen.core.structure")
    c = s.Structure.from_file(cif_path)
    num_atoms = size(c)[1]

    # list of atom symbols
    atom_ids = [site_element(s) for s in c]

    if use_voronoi
        weight_mat = build_graph_voronoi(c, num_atoms)
    else
        weight_mat = build_graph_cutoff(c, radius, max_num_nbr, dist_decay_func)
    end

    # normalize weights
    if normalize
        weight_mat = weight_mat ./ maximum(weight_mat)
    end

    return (weight_mat, atom_ids)
end

"""
Build graph using neighbors from faces of Voronoi polyedra and weights from areas. Based on the approach from https://github.com/ulissigroup/uncertainty_benchmarking/blob/aabb407807e35b5fd6ad06b14b440609ae09e6ef/BNN/data_pyro.py#L268
"""
function build_graph_voronoi(struc)
    num_atoms = size(struc)[1]
    sa = pyimport("pymatgen.analysis.structure_analyzer")
    vc = sa.VoronoiConnectivity(struc)
    conn = vc.connectivity_array
    weight_mat = zeros(Float32, num_atoms, num_atoms)
    # loop over central atoms
    for atom_ind in 1:size(conn)[1]
        # loop over neighbor atoms
        for nb_ind in 1:size(conn)[2]
            # loop over each possible PBC image for chosen image
            for image_ind in 1:size(conn)[3]
                # only add as neighbor if atom is not current center one AND there is connectivity to image
                if (atom_ind != image_ind) & (conn[atom_ind, nb_ind, image_ind] != 0)
                    weight_mat[atom_ind, nb_ind] += conn[atom_ind, nb_ind, image_ind]/maximum(conn[atom_ind, :, :])
                end
            end
        end
    end

    # average across diagonal (because neighborness isn't strictly symmetric in the way we're defining it here)
    weight_mat = 0.5.* (weight_mat .+ weight_mat')
end

"""
Build graph using neighbor number cutoff method adapted from original CGCNN. Note that `max_num_nbr` is a "soft" max, in that if there are more of the same distance as the last, all of those will be added.

# TODO
- option to cut off by nearest, next-nearest, etc. by DISTANCE rather than NUMBER of neighbors
"""
function build_graph_cutoff(struc; radius=8.0, max_num_nbr=12, dist_decay_func=inverse_square)
    #= find neighbors, requires a cutoff radius
    returns a NxM Array of PyObject PeriodicSite
    ...except when it returns a list of N of lists of length M...
    each PeriodicSite is (site, distance, index, image)
    N is num sites in crystal, M is num neighbors
    =#
    all_nbrs = c.get_all_neighbors(radius, include_index=true)

    # sort by distance
    # returns list of length N of lists of length M
    if length(size(all_nbrs)) == 2
        all_nbrs = [sort(all_nbrs[i,:], lt=(x,y)->isless(site_distance(x), site_distance(y))) for i in 1:num_atoms]
    elseif length(size(all_nbrs)) == 1
        all_nbrs = [sort(all_nbrs[i][:], lt=(x,y)->isless(site_distance(x), site_distance(y))) for i in 1:num_atoms]
    end

    # iterate through each list of neighbors (corresponding to neighbors of a given atom) to find bonds (eventually, graph edges)
    weight_mat = zeros(Float32, num_atoms, num_atoms)
    for atom_ind in 1:num_atoms
        this_atom = get(c, atom_ind-1)
        atom_nbs = all_nbrs[atom_ind]
        # iterate over each neighbor...
        for nb_num in 1:size(all_nbrs[atom_ind])[1]
            nb = atom_nbs[nb_num]
            global nb_ind = site_index(nb)
            # if we're under the max, add it for sure
            if nb_num < max_num_nbr
                weight_mat[atom_ind, nb_ind] = weight_mat[atom_ind, nb_ind] + dist_decay_func(site_distance(nb))
            # if we're at/above the max, add if distance is the same
            else
                # check we're not on the last one
                if nb_ind < size(atom_nbs)[1] - 1
                    next_nb = atom_nbs[nb_ind + 1]
                    # add another bond if it's the exact same distance to the next neighbor in the list
                    if are_equidistant(nb, next_nb)
                        weight_mat[atom_ind, nb_ind] = weight_mat[atom_ind, nb_ind] + dist_decay_func(site_distance(nb))
                    end
                end
            end
        end
    end

    # average across diagonal (because neighborness isn't strictly symmetric in the way we're defining it here)
    weight_mat = 0.5.* (weight_mat .+ weight_mat')
end
