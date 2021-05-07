using LightGraphs
using SimpleWeightedGraphs
using LinearAlgebra
using GraphPlot
using Colors

# TO CONSIDER: store ref to featurization rather than the thing itself? Does this matter for any performance we care about?
"""
    AtomGraph

A type representing an atomic structure as a graph (`gr`).

# Fields
- `graph::SimpleWeightedGraph{Int32,Float32}`: the graph representing the structure. See
  [`build_graph`](@ref) for more on generating the weights.
- `elements::Vector{String}`: list of elemental symbols corresponding to each node of the
  graph
- `laplacian::Matrix{Float32}`: Normalized graph Laplacian matrix, stored to speed up
  convolution operations by avoiding recomputing it every pass.
- `features::Matrix{Float32}`: Feature matrix of size (# features, # nodes). AtomGraph can
  be initialized without defining this field, but if it is defined, the subsequent field must be also.
- `featurization::AbstractFeaturization`: Featurization scheme specification to maintain 
  "decodability" of features.
- `id::String`: Optional, an identifier, e.g. to correspond with tags/labels of an imported
  dataset.
"""
mutable struct AtomGraph <: AbstractAtoms
    graph::SimpleWeightedGraph{<: Integer, <: Real}
    elements::Vector{String}
    laplacian::Matrix{<: Real} # wanted to use LightGraphs.LinAlg.NormalizedGraphLaplacian but seems this doesn't support weighted graphs?
    atom_features::Union{Matrix{<: Real},Nothing} # if we add edge features this type will have to relax
    featurization::Union{AbstractFeaturization,Nothing}
    id::String # or maybe we let it be a number too?
end


# first, the basic constructor
"""
    AtomGraph(gr, elements, features, featurization, id="")
    AtomGraph(gr, elements, id="")
    AtomGraph(adj, elements, features, featurization, id="")
    AtomGraph(adj, elements, id="")

Construct an AtomGraph object, either directly from a SimpleWeightedGraph `gr` or from an
adjacency matrix `adj`, along with, at minimum, the list of elemental symbols `elements`
representing each node. Note that the object can be initialized without features, but if
features are provided, so too must be the featurization scheme, in order to maintain
"decodability" of features.
"""
function AtomGraph(
    gr::SimpleWeightedGraph{A, B},
    elements::Vector{String},
    features::Matrix{<: Real},
    featurization::AbstractFeaturization,
    id = "",
) where {B <: Real, A <: Integer}
    # check that elements is the right length
    num_atoms = size(gr)[1]
    @assert length(elements) == num_atoms "Element list length doesn't match graph size!"

    # check that features is the right dimensions (# features x # nodes)
    expected_feature_length = sum(f.num_bins for f in featurization)
    @assert size(features) == (expected_feature_length, num_atoms) "Feature matrix is of wrong dimension! It should be of size (# features, # nodes)"

    # if all these are good, calculate laplacian and build the thing
    laplacian = normalized_laplacian(gr)
    AtomGraph(gr, elements, laplacian, features, featurization, id)
end


# one without features or featurization initialized yet
function AtomGraph(
    gr::SimpleWeightedGraph{A, B},
    elements::Vector{String},
    id = ""
) where {B <: Real,A <: Integer}
    # check that elements is the right length
    num_atoms = size(gr)[1]
    @assert length(elements) == num_atoms "Element list length doesn't match graph size!"

    laplacian = B.(normalized_laplacian(gr))
    AtomGraph(gr, elements, laplacian, nothing, nothing, id)
end


# initialize directly from adjacency matrix
AtomGraph(
    adj::Array{R},
    elements::Vector{String},
    features::Matrix{R},
    featurization::AbstractFeaturization,
    id = ""
) where {R <: Real} =
    AtomGraph(SimpleWeightedGraph{R}(adj), elements, features, featurization, id)


AtomGraph(adj::Array{R}, elements::Vector{String}, id = "") where {R <: Real} =
    AtomGraph(SimpleWeightedGraph{R}(adj), elements, id)


# pretty printing, short version
function Base.show(io::IO, ag::AtomGraph)
    st = "AtomGraph $(ag.id) with $(nv(ag.graph)) nodes, $(ne(ag.graph)) edges"
    if !isnothing(ag.featurization)
        st = string(st, ", feature vector length $(size(ag.atom_features)[1])")
    end
    print(io, st)
end


# pretty printing, long version
function Base.show(io::IO, ::MIME"text/plain", ag::AtomGraph)
    st = "AtomGraph $(ag.id) with $(nv(ag.graph)) nodes, $(ne(ag.graph)) edges\n   atoms: $(ag.elements)\n   feature vector length: "
    if isnothing(ag.featurization)
        st = string(st, "uninitialized\n   encoded features: uninitialized")
    else
        st =
            string(st, "$(size(ag.atom_feats)[1])\n   featurization: ", ag.featurization)
    end
    print(io, st)

end


"""
    normalized_laplacian(graph)

Compute the normalized graph Laplacian matrix of the input graph, defined as

``I - D^{-1/2} A D^{-1/2}``

where ``A`` is the adjacency matrix and ``D`` is the degree matrix.
"""
function normalized_laplacian(g::G) where {G <: LightGraphs.AbstractGraph}
    a = adjacency_matrix(g)
    d = vec(sum(a, dims = 1))
    inv_sqrt_d = diagm(0 => d .^ (-0.5f0))
    laplacian = Float32.(I - inv_sqrt_d * a * inv_sqrt_d)
    !any(isnan, laplacian) || throw(
        ArgumentError(
            "NaN values in graph Laplacian! This is most likely due to atomic separations larger than the specified cutoff distance leading to block zeros in the adjacency matrix...try increasing the cutoff distance or inspecting your structure to ensure the file is correct.",
        ),
    )
    return laplacian
end


normalized_laplacian(ag::AtomGraph) = ag.laplacian


# now visualization stuff...

"Get a list of colors to use for graph visualization."
function graph_colors(atno_list, seed_color = colorant"cyan4")
    atom_types = unique(atno_list)
    atom_type_inds = Dict(atom_types[i] => i for i = 1:length(atom_types))
    color_inds = [atom_type_inds[i] for i in atno_list]
    colors = distinguishable_colors(length(atom_types), seed_color)
    return colors[color_inds]
end


# helper fcn for sorting because edge ordering isn't preserved when converting to SimpleGraph
function lt_edge(
    e1::SimpleWeightedGraphs.SimpleWeightedEdge{<: Integer, <: Real},
    e2::SimpleWeightedGraphs.SimpleWeightedEdge{<: Integer, <: Real},
)
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
    edges_sorted = sort([e for e in edges(ag.graph)], lt = lt_edge)
    for e in edges_sorted
        append!(edgewidths, e.weight)
    end
    return edgewidths
end


"Visualize a given graph."
function visualize(ag::AtomGraph)
    # gplot doesn't work on weighted graphs
    sg = SimpleGraph(adjacency_matrix(ag))
    plt = gplot(
        sg,
        nodefillc = graph_colors(ag.elements),
        nodelabel = ag.elements,
        edgelinewidth = graph_edgewidths(ag),
    )
    display(plt)
end
