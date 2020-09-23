using LightGraphs; const lg=LightGraphs
using SimpleWeightedGraphs
using LinearAlgebra
using GraphPlot
using Colors
using JSON

# Type to store atomic graphs
# TO CONSIDER: store ref to featurization rather than the thing itself? Does this matter for any performance we care about?
# TO CONSIDER: store an ID of some kind? (e.g. mp-123 or whatever)
mutable struct AtomGraph <: lg.AbstractGraph{Float32}
    graph::SimpleWeightedGraph{Int32,Float32}
    elements::Vector{String} # list of elemental symbols corresponding to each node
    lapl::Matrix{Float32} # graph laplacian (normalized)
    features::Matrix{Float32} # feature matrix (size (# features, # nodes))
    featurization::Vector{AtomFeat} # featurization scheme in the form of a list of AtomFeat objects
end

# basic constructor
function AtomGraph(gr::SimpleWeightedGraph{Int32,Float32}, el_list::Vector{String}, features::Matrix{Float32}, featurization::Vector{AtomFeat})
    # check that el_list is the right length
    num_atoms = size(gr)[1]
    @assert length(el_list)==num_atoms "Element list length doesn't match graph size!"

    # check that features is the right dimensions (# features x # nodes)
    expected_feature_length = sum([f.num_bins for f in featurization])
    @assert size(features) == (expected_feature_length, num_atoms) "Feature matrix is of wrong dimension! It should be of size (# features, # nodes)"

    # if all these are good, calculate laplacian and build the thing
    lapl = normalized_laplacian(gr)
    AtomGraph(gr, el_list, lapl, features, featurization)
end

# one without features or featurization initialized yet
function AtomGraph(gr::SimpleWeightedGraph{Int32,Float32}, el_list::Vector{String})
    # check that el_list is the right length
    num_atoms = size(gr)[1]
    @assert length(el_list)==num_atoms "Element list length doesn't match graph size!"

    lapl = Float32.(normalized_laplacian(gr))
    AtomGraph(gr, el_list, lapl, zeros(Float32, 1, num_atoms), AtomFeat[])
end

# TODO, maybe: constructor where you give adjacency matrix and it builds the graph for you also

JSON.lower(ag::AtomGraph) = Dict("graph"=>ag.graph.weights, "elements"=>ag.elements, "lapl"=>ag.lapl, "features"=>ag.features, "featurization"=>JSON.lower.(ag.featurization))

function save_json(ag::AtomGraph, fpath::String)
    open(fpath, "w") do f
        write(f, json(ag))
    end
end

function AtomGraph(d::Dict{String,Any})
    graph = SimpleWeightedGraph{Int32,Float32}(Float32.(hcat([v for v in d["graph"]]...)))
    elements = d["elements"]
    lapl = Float32.(hcat([v for v in d["lapl"]]...))
    features = Float32.(hcat([v for v in d["features"]]...))
    featurization = [AtomFeat(fd) for fd in d["featurization"]]
end

AtomGraph(json_path::String) = AtomGraph(JSON.parsefile(json_path))

# pretty printing, short version
function Base.show(io::IO, g::AtomGraph)
    st = "AtomGraph with $(nv(g)) nodes, $(ne(g)) edges"
    if length(g.featurization)!=0
        st = string(st, ", feature vector length $(size(g.features)[1])")
    end
    print(io, st)
end

# pretty printing, long version
function Base.show(io::IO, ::MIME"text/plain", g::AtomGraph)
    st = "AtomGraph with $(nv(g)) nodes, $(ne(g)) edges\n   atoms: $(g.elements)\n   feature vector length: "
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

# TODO: maybe some subgraph stuff for cutting up graphs later, will need to make sure to account for feature matrix and element list properly as well as recompute laplacian, would be cool to be able to e.g. filter on particular features and only pull matching nodes

# this cribbed from GeometricFlux
function normalized_laplacian(g::G) where G<:lg.AbstractGraph
    a = adjacency_matrix(g)
    d = vec(sum(a, dims=1))
    inv_sqrt_d = diagm(0=>d.^(-0.5))
    Float32.(I - inv_sqrt_d * a * inv_sqrt_d)
end

normalized_laplacian(g::AtomGraph) = g.lapl

# function to add node features if only the "bare" graph has been initialized, note that featurization scheme must be specified!
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
    feature_mat, featurization = make_feature_vectors(featurization)
    add_features!(g, feature_mat, featurization)
end

# is there a clever way to roll this into the previous one since the syntax is identical?
function add_features!(g::AtomGraph, feature_names::Vector{Symbol})
    feature_mat, featurization = make_feature_vectors(feature_names)
    add_features!(g, feature_mat, featurization)
end

function add_features!(g::AtomGraph, feature_names::Vector{Symbol}, nbins::Vector{<:Integer}, logspaced=false)
    feature_mat, featurization = make_feature_vectors(feature_names, nbins=nbins, logspaced=logspaced)
    add_features!(g, feature_mat, featurization)
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

"Compute edge widths (proportional to weights on graph) for graph visualization."
function graph_edgewidths(g, weight_mat)
    edgewidths = []
    # should be able to do this as
    for e in edges(g)
        append!(edgewidths, weight_mat[e.src, e.dst])
    end
    return edgewidths
end

"Visualize a given graph."
function visualize_graph(g, element_list)
    # gplot doesn't work on weighted graphs
    sg = SimpleGraph(adjacency_matrix(g))
    plt = gplot(sg, nodefillc=graph_colors(element_list), nodelabel=element_list, edgelinewidth=graph_edgewidths(sg, g.weights))
    display(plt)
end