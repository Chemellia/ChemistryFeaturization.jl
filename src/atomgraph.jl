using LightGraphs; const lg=LightGraphs
using SimpleWeightedGraphs
using LinearAlgebra
using GraphPlot
using Colors
using JSON

# TO CONSIDER: store ref to featurization rather than the thing itself? Does this matter for any performance we care about?
"""
    AtomGraph

A type representing an atomic structure as a graph (`gr`).

# Fields
- `graph::SimpleWeightedGraph{Int32,Float32}`: the graph representing the structure. See
  [`build_graph`](@ref) for more on generating the weights.
- `elements::Vector{String}`: list of elemental symbols corresponding to each node of the
  graph
- `lapl::Matrix{Float32}`: Normalized graph Laplacian matrix, stored to speed up
  convolution operations by avoiding recomputing it every pass.
- `features::Matrix{Float32}`: Feature matrix of size (# features, # nodes). AtomGraph can
  be initialized without defining this field, but if it is defined, the subsequent field must be also.
- `featurization::Vector{AtomFeat}`: Featurization scheme specification to maintain 
  "decodability" of features.
- `id::String`: Optional, an identifier, e.g. to correspond with tags/labels of an imported
  dataset.
"""
mutable struct AtomGraph <: lg.AbstractGraph{Float32}
    graph::SimpleWeightedGraph{Int32,Float32}
    elements::Vector{String}
    lapl::Matrix{Float32}
    features::Matrix{Float32}
    featurization::Vector{AtomFeat}
    id::String
end

# first, the basic constructor
"""
    AtomGraph(gr, el_list, features, featurization, id="")
    AtomGraph(gr, el_list, id="")
    AtomGraph(adj, el_list, features, featurization, id="")
    AtomGraph(adj, el_list, id="")

Construct an AtomGraph object, either directly from a SimpleWeightedGraph `gr` or from an
adjacency matrix `adj`, along with, at minimum, the list of elemental symbols `el_list`
representing each node. Note that the object can be initialized without features, but if
features are provided, so too must be the featurization scheme, in order to maintain
"decodability" of features.
"""
function AtomGraph(gr::SimpleWeightedGraph{Int32,Float32}, el_list::Vector{String}, features::Matrix{Float32}, featurization::Vector{AtomFeat}, id="")
    # check that el_list is the right length
    num_atoms = size(gr)[1]
    @assert length(el_list)==num_atoms "Element list length doesn't match graph size!"

    # check that features is the right dimensions (# features x # nodes)
    expected_feature_length = sum(f.num_bins for f in featurization)
    @assert size(features) == (expected_feature_length, num_atoms) "Feature matrix is of wrong dimension! It should be of size (# features, # nodes)"

    # if all these are good, calculate laplacian and build the thing
    lapl = normalized_laplacian(gr)
    AtomGraph(gr, el_list, lapl, features, featurization, id)
end

# one without features or featurization initialized yet
function AtomGraph(gr::SimpleWeightedGraph{Int32,Float32}, el_list::Vector{String}, id="")
    # check that el_list is the right length
    num_atoms = size(gr)[1]
    @assert length(el_list)==num_atoms "Element list length doesn't match graph size!"

    lapl = Float32.(normalized_laplacian(gr))
    AtomGraph(gr, el_list, lapl, zeros(Float32, 1, num_atoms), AtomFeat[], id)
end

# initialize directly from adjacency matrix
AtomGraph(adj::Array{Float32}, el_list::Vector{String}, features::Matrix{Float32}, featurization::Vector{AtomFeat}, id="") = AtomGraph(SimpleWeightedGraph{Int32}(adj), el_list, features, featurization, id)
AtomGraph(adj::Array{Float32}, el_list::Vector{String}, id="") = AtomGraph(SimpleWeightedGraph{Int32}(adj), el_list, id)

# pretty printing, short version
function Base.show(io::IO, g::AtomGraph)
    st = "AtomGraph $(g.id) with $(nv(g)) nodes, $(ne(g)) edges"
    if length(g.featurization)!=0
        st = string(st, ", feature vector length $(size(g.features)[1])")
    end
    print(io, st)
end

# pretty printing, long version
function Base.show(io::IO, ::MIME"text/plain", g::AtomGraph)
    st = "AtomGraph $(g.id) with $(nv(g)) nodes, $(ne(g)) edges\n   atoms: $(g.elements)\n   feature vector length: "
    if length(g.featurization)==0
        st = string(st, "uninitialized\n   encoded features: uninitialized")
    else
        st = string(st, "$(size(g.features)[1])\n   encoded features: ",  [string(f.name, ", ") for f in g.featurization]...)[1:end-2]
    end
    print(io, st)
end

# now define some functions...

# first, the LightGraphs required ones (easy, just run them on the graph contents)
lg.edges(g::AtomGraph) = lg.edges(g.graph)
lg.eltype(g::AtomGraph) = lg.eltype(g.graph)
lg.edgetype(g::AtomGraph) = lg.edgetype(g.graph)
lg.ne(g::AtomGraph) = lg.ne(g.graph)
lg.nv(g::AtomGraph) = lg.nv(g.graph)
lg.weights(g::AtomGraph) = lg.weights(g.graph)
lg.is_directed(g::AtomGraph) = lg.is_directed(g.graph)

lg.outneighbors(g::AtomGraph, node) = lg.outneighbors(g.graph, node)
lg.inneighbors(g::AtomGraph, node) = lg.inneighbors(g.graph, node)
lg.has_vertex(g::AtomGraph, v::Integer) = lg.has_vertex(g.graph, v)
lg.has_edge(g::AtomGraph, i, j) = lg.has_edge(g.graph, i, j)

lg.zero(AtomGraph) = AtomGraph(zero(SimpleWeightedGraph{Int32,Float32}), String[])

"""
    normalized_laplacian(graph)

Compute the normalized graph Laplacian matrix of the input graph, defined as

``I - D^{-1/2} A D^{-1/2}``

where ``A`` is the adjacency matrix and ``D`` is the degree matrix.
"""
function normalized_laplacian(g::G) where G<:lg.AbstractGraph
    a = adjacency_matrix(g)
    d = vec(sum(a, dims=1))
    inv_sqrt_d = diagm(0=>d.^(-0.5f0))
    lapl = Float32.(I - inv_sqrt_d * a * inv_sqrt_d)
    !any(isnan, lapl) || throw(ArgumentError("NaN values in graph Laplacian! This is most likely due to atomic separations larger than the specified cutoff distance leading to block zeros in the adjacency matrix...try increasing the cutoff distance or inspecting your structure to ensure the file is correct."))
    return lapl
end

normalized_laplacian(g::AtomGraph) = g.lapl

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

# now visualization stuff...

"Get a list of colors to use for graph visualization."
function graph_colors(atno_list, seed_color=colorant"cyan4")
    atom_types = unique(atno_list)
    atom_type_inds = Dict(atom_types[i]=>i for i in 1:length(atom_types))
    color_inds = [atom_type_inds[i] for i in atno_list]
    colors = distinguishable_colors(length(atom_types), seed_color)
    return colors[color_inds]
end

# helper fcn for sorting because edge ordering isn't preserved when converting to SimpleGraph
function lt_edge(e1::SimpleWeightedGraphs.SimpleWeightedEdge{Int32,Float32}, e2::SimpleWeightedGraphs.SimpleWeightedEdge{Int32,Float32})
    if e1.src < e2.src
        return true
    elseif e1.dst < e2.dst
        return true
    else
        return false
    end
end

"Compute edge widths (proportional to weights on graph) for graph visualization."
function graph_edgewidths(ag::AtomGraph)
    edgewidths = []
    edges_sorted = sort([e for e in edges(ag)], lt=lt_edge)
    for e in edges_sorted
        append!(edgewidths, e.weight)
    end
    return edgewidths
end

"Visualize a given graph."
function visualize_graph(ag::AtomGraph)
    # gplot doesn't work on weighted graphs
    sg = SimpleGraph(adjacency_matrix(ag))
    plt = gplot(sg, nodefillc=graph_colors(ag.elements), nodelabel=ag.elements, edgelinewidth=graph_edgewidths(ag))
    display(plt)
end
