# the rest of this will either be unnecessary and vanish, or change and move to utils for automatically generating encode_f and decode_f stuff

"""
    add_features!(ag, features, featurization)
    add_features!(ag, atom_feature_vecs, featurization)
    add_features!(ag, featurization)
    add_features!(ag, feature_names; nbins, logspaced=false)

Add atomic features to an existing AtomGraph object. Can be done by (in decreasing order of
speed since more things have to be calculated; these bullet points correspond to the four
call signature options above):
- providing feature matrix `features` (size # features x # nodes) directly
- providing a dictionary of feature vectors (keys should be elemental symbols) and matrix
  will be built
- providing only a featurization scheme (list of AtomFeat objects), and vectors for each 
  element will be built and combined into a feature matrix
- Providing a list of feature names (and optionally specify binning and spacing for 
  continuous numerical features) and built-in data will be used to build featurization
  scheme, vectors, and matrix.

In every case, a featurization scheme (or ingredients to build it) must be provided to 
ensure decodability of features.

See also: [`add_features_batch!`](@ref)
"""
function add_features!(g::AtomGraph, features::Matrix{Float32}, featurization::Vector{AtomFeat})
    num_atoms = nv(g)

    # check that features is the right dimensions (# features x # nodes)
    expected_feature_length = sum([f.num_bins for f in featurization])
    @assert size(features) == (expected_feature_length, num_atoms) "Feature matrix is of wrong dimension! It should be of size (# features, # nodes)"

    # okay now we can set the features
    g.features = features
    g.featurization = featurization
end

# alternate version where it builds the features too, you have to pass in the results of the make_feature_vectors function
function add_features!(g::AtomGraph, atom_feature_vecs::Dict{String, Vector{Float32}}, featurization::Vector{AtomFeat})
    @assert Set(String.(g.elements)) <= Set(keys(atom_feature_vecs)) "Some atoms in your graph do not have corresponding feature vectors! This could be because some features you requested had missing values for these atoms."
    feature_mat = Float32.(hcat([atom_feature_vecs[e] for e in g.elements]...))
    add_features!(g, feature_mat, featurization)
end

# and finally the ones where it makes the feature vectors too...(defined for both signatures of the make_feature_vectors function just for completeness)
function add_features!(g::AtomGraph, featurization::Vector{AtomFeat})
    feature_vecs, featurization = make_feature_vectors(featurization)
    add_features!(g, feature_vecs, featurization)
end

#function add_features!(g::AtomGraph, feature_names::Vector{Symbol})
#    feature_vecs, featurization = make_feature_vectors(feature_names)
#    add_features!(g, feature_vecs, featurization)
#end

function add_features!(g::AtomGraph, feature_names::Vector{Symbol}; nbins::Vector{<:Integer}, logspaced=false)
    feature_vecs, featurization = make_feature_vectors(feature_names, nbins=nbins, logspaced=logspaced)
    add_features!(g, feature_vecs, featurization)
end

"""
    add_features_batch!(ags, atom_feature_vecs, featurization)
    add_features_batch!(ags, featurization)
    add_features_batch!(ags, feature_names; nbins, logspaced=false)

Add features to a list of AtomGraph objects. See [`add_features!`](@ref).
"""
function add_features_batch!(gs::Array{AtomGraph}, atom_feature_vecs::Dict{String, Vector{Float32}}, featurization::Vector{AtomFeat})
    for g in gs
        add_features!(g, atom_feature_vecs, featurization)
    end
end

function add_features_batch!(gs::Array{AtomGraph}, featurization::Vector{AtomFeat})
    feature_vecs, featurization = make_feature_vectors(featurization)
    add_features_batch!(gs, feature_vecs, featurization)
end

function add_features_batch!(gs::Array{AtomGraph}, feature_names::Vector{Symbol}; nbins::Vector{<:Integer}, logspaced=false)
    feature_vecs, featurization = make_feature_vectors(feature_names, nbins=nbins, logspaced=logspaced)
    add_features_batch!(gs, feature_vecs, featurization)
end
