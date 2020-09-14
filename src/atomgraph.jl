using LightGraphs; const lg=LightGraphs
using SimpleWeightedGraphs
using LinearAlgebra

# Type to store atomic graphs
# TO CONSIDER: store ref to featurization rather than the thing itself? Does this matter?
mutable struct AtomGraph{T <: AbstractSimpleWeightedGraph{Int64, Float32}} <: lg.AbstractGraph{Float32}
    graph::T
    elements::Vector{String}
    lapl::Matrix{Float32}
    features::Matrix{Float32}
    featurization::Vector{AtomFeat}
end

# basic constructor
function AtomGraph(gr::G, el_list::Vector{String}; features::Matrix{Float32}, featurization::Vector{AtomFeat}) where G <: lg.AbstractGraph
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
function AtomGraph(gr::G, el_list::Vector{String}) where G <: lg.AbstractGraph
    # check that el_list is the right length
    num_atoms = size(gr)[1]
    @assert length(el_list)==num_atoms "Element list length doesn't match graph size!"

    lapl = Float32.(normalized_laplacian(gr))
    AtomGraph(gr, el_list, lapl, zeros(Float32, 1, num_atoms), AtomFeat[])
end

# TODO, maybe: constructor where you give adjacency matrix and it builds the graph for you also

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

# TODO: implement zero to return one with zero(G) and use constructors from above to populate other stuff, this is the last one for all of the LightGraphs functions to "just work"

# this cribbed from GeometricFlux
function normalized_laplacian(g::G) where G<:lg.AbstractGraph
    a = adjacency_matrix(g)
    d = vec(sum(a, dims=1))
    inv_sqrt_d = diagm(0=>d.^(-0.5))
    I - inv_sqrt_d * a * inv_sqrt_d
end

# function to add node features if only the "bare" graph has been initialized
# featurization scheme must be specified!
function add_features!(g::AtomGraph, features::Matrix{Float32}, featurization::Vector{AtomFeat})
    num_atoms = nv(g)

    # check that features is the right dimensions (# features x # nodes)
    expected_feature_length = sum([f.num_bins for f in featurization])
    @assert size(features) == (expected_feature_length, num_atoms) "Feature matrix is of wrong dimension! It should be of size (# features, # nodes)"

    # okay now we can set the features
    g.features = features
    g.featurization = featurization
end

# TODO: pretty printing: https://docs.julialang.org/en/v1/manual/types/#man-custom-pretty-printing

# TODO: visualization: pull in stuff from graph_vis.jl