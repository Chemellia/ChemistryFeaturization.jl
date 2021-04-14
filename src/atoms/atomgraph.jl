using LightGraphs
using SimpleWeightedGraphs
using LinearAlgebra
using GraphPlot
using Colors
using ..ChemistryFeaturization: GraphNodeFeaturization

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
mutable struct AtomGraph{T<:Real} <: AbstractAtoms
    graph::SimpleWeightedGraph{<:Unsigned,T}
    elements::Vector{String}
    lapl::Matrix{T} # wanted to use LightGraphs.LinAlg.NormalizedGraphLaplacian but seems this doesn't support weighted graphs?
    atom_feats::Matrix{<:Real} # if we add edge features this type will have to relax
    featurization::GraphNodeFeaturization
    id::String # or maybe we let it be a number too?
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
function AtomGraph(gr::SimpleWeightedGraph{<:Unsigned,<:Real}, el_list::Vector{String}, features::Matrix{<:Real}, featurization::GraphNodeFeaturization, id="")
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
function AtomGraph(gr::SimpleWeightedGraph{A<:Unsigned,B<:Real}, el_list::Vector{String}, id="")
    # check that el_list is the right length
    num_atoms = size(gr)[1]
    @assert length(el_list)==num_atoms "Element list length doesn't match graph size!"

    lapl = B.(normalized_laplacian(gr))
    AtomGraph(gr, el_list, lapl, zeros(B, 1, num_atoms), AtomFeat[], id)
end

# initialize directly from adjacency matrix, defaults to 32-bit integers because unlikely to need more nodes than that...
AtomGraph(adj::Array{<:Real}, el_list::Vector{String}, features::Matrix{<:Real}, featurization::GraphNodeFeaturization, id=""; U=UInt32) = AtomGraph(SimpleWeightedGraph{T}(adj), el_list, features, featurization, id)
AtomGraph(adj::Array{<:Real}, el_list::Vector{String}, id=""; U=UInt32) = AtomGraph(SimpleWeightedGraph{T}(adj), el_list, id)

# pretty printing, short version
function Base.show(io::IO, ag::AtomGraph)
    st = "AtomGraph $(ag.id) with $(nv(ag.graph)) nodes, $(ne(ag.graph)) edges"
    if length(g.featurization)!=0
        st = string(st, ", feature vector length $(size(ag.atom_feats)[1])")
    end
    print(io, st)
end

# pretty printing, long version
function Base.show(io::IO, ::MIME"text/plain", ag::AtomGraph)
    st = "AtomGraph $(ag.id) with $(nv(ag.graph)) nodes, $(ne(ag.graph)) edges\n   atoms: $(ag.elements)\n   feature vector length: "
    if length(g.featurization)==0
        st = string(st, "uninitialized\n   encoded features: uninitialized")
    else
        st = string(st, "$(size(ag.features)[1])\n   encoded features: ",  [string(f.name, ", ") for f in g.featurization]...)[1:end-2]
    end
    print(io, st)
end

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

normalized_laplacian(ag::AtomGraph) = ag.lapl

# maybe some cutesy stuff like dispatching things like length as length(elements)

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
function lt_edge(e1::SimpleWeightedGraphs.SimpleWeightedEdge{<:Unsigned,<:Real}, e2::SimpleWeightedGraphs.SimpleWeightedEdge{<:Unsigned,<:Real})
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
    edges_sorted = sort([e for e in edges(ag.graph)], lt=lt_edge)
    for e in edges_sorted
        append!(edgewidths, e.weight)
    end
    return edgewidths
end

"Visualize a given graph."
function visualize(ag::AtomGraph)
    # gplot doesn't work on weighted graphs
    sg = SimpleGraph(adjacency_matrix(ag))
    plt = gplot(sg, nodefillc=graph_colors(ag.elements), nodelabel=ag.elements, edgelinewidth=graph_edgewidths(ag))
    display(plt)
end