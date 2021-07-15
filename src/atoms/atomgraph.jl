using LightGraphs
using SimpleWeightedGraphs
using LinearAlgebra
using GraphPlot
using Colors
using Serialization
using ..ChemistryFeaturization.Utils.GraphBuilding

using ..ChemistryFeaturization.AbstractType: AbstractAtoms

# TO CONSIDER: store ref to featurization rather than the thing itself? Does this matter for any performance we care about?
"""
    AtomGraph

A type representing an atomic structure as a graph (`gr`).

# Fields
- `graph::SimpleWeightedGraph{<:Integer,<:Real}`: the graph representing the structure. See
  [`build_graph`](@ref) for more on generating the weights.
- `elements::Vector{String}`: list of elemental symbols corresponding to each node of the
  graph
- `laplacian::Matrix{<:Real}`: Normalized graph Laplacian matrix, stored to speed up
  convolution operations by avoiding recomputing it every pass.
- `id::String`: Optional, an identifier, e.g. to correspond with tags/labels of an imported
  dataset.
"""
mutable struct AtomGraph <: AbstractAtoms
    graph::SimpleWeightedGraph{<:Integer,<:Real}
    elements::Vector{String}
    laplacian::Matrix{<:Real} # wanted to use LightGraphs.LinAlg.NormalizedGraphLaplacian but seems this doesn't support weighted graphs?
    id::String # or maybe we let it be a number too?
end

# one without features or featurization initialized yet
function AtomGraph(
    graph::SimpleWeightedGraph{A,B},
    elements::Vector{String},
    id = "",
) where {B<:Real,A<:Integer}
    # check that elements is the right length
    num_atoms = size(graph)[1]
    @assert length(elements) == num_atoms "Element list length doesn't match graph size!"

    # this was previously B.(normalized_laplacian(graph)) - won't that potentially give rise to compatibility issues if B is a custom type?
    laplacian = normalized_laplacian(graph)
    AtomGraph(graph, elements, laplacian, id)
end


# initialize directly from adjacency matrix
AtomGraph(adj::Array{R}, elements::Vector{String}, id = "") where {R<:Real} =
    AtomGraph(SimpleWeightedGraph(adj), elements, id)

"""
    AtomGraph(input_file_path, id = splitext(input_file_path)[begin]; output_file_path = nothing, featurization = nothing, overwrite_file = false, use_voronoi = false, cutoff_radius = 8.0, max_num_nbr = 12, dist_decay_func = inverse_square, normalize_weights = true)

Construct an AtomGraph object from a structure file.

# Required Arguments
- `input_file_path::String`: path to file containing structure (must be readable by ASE.io.read)

# Optional Arguments
- `id::String`: ID associated with structure (e.g. identifier from online database). Defaults to name of input file if undefined.
- `output_file_path = nothing`: If provided, structure will be serialized to file at this location
- `overwrite_file::Bool = false`: whether to overwrite an existing file at `output_file_path`
- `use_voronoi::Bool = false`: Whether to build neighbor lists using Voronoi decompositions
- `cutoff_radius::Real = 8.0`: If not using Voronoi neighbor lists, longest allowable distance to a neighbor, in Angstroms
- `max_num_nbr::Integer = 12`: If not using Voronoi neighbor lists, largest allowable number of neighbors
- `dist_decay_func = inverse_square`: Function by which to assign edge weights according to distance between neighbors
- `normalize_weights::Bool = true`: Whether to normalize weights such that the largest is 1.0

# Note
`max_num_nbr` is a "soft" limit – if multiple neighbors are at the same distance, the full neighbor list may be longer.
"""
function AtomGraph(
    input_file_path::String,
    id::String = splitext(input_file_path)[begin];
    output_file_path::Union{String,Nothing} = nothing,
    overwrite_file::Bool = false,
    use_voronoi::Bool = false,
    cutoff_radius::Real = 8.0,
    max_num_nbr::Integer = 12,
    dist_decay_func::Function = inverse_square,
    normalize_weights::Bool = true,
)

    local ag

    if !isfile(input_file_path)
        @warn "$input_file_path does not exist. Cannot build graph from a non-existent file."
        return missing
    end

    if splitext(input_file_path)[end] == ".jls" # deserialize
        ag = deserialize(input_file_path)
        ag.id = id

    else # try actually building the graph
        try
            adj_mat, elements = build_graph(
                input_file_path,
                use_voronoi = use_voronoi,
                cutoff_radius = cutoff_radius,
                max_num_nbr = max_num_nbr,
                dist_decay_func = dist_decay_func,
                normalize_weights = normalize_weights,
            )
            ag = AtomGraph(adj_mat, elements, id)
        catch
            @warn "Unable to build graph for $input_file_path"
            return missing
        end
    end

    to_serialize = !isnothing(output_file_path)
    if to_serialize
        if isfile(output_file_path) && !(overwrite_file)
            @info "Output file already exists, and `overwrite_file` is set to false.\nIf you want to overwrite the existing graph, set `overwrite=true`, or remove the existing file and retry."
        else
            serialize(output_file_path, ag)
        end
    end

    return ag
end

# pretty printing, short version
function Base.show(io::IO, ag::AtomGraph)
    st = "AtomGraph $(ag.id) with $(nv(ag.graph)) nodes, $(ne(ag.graph)) edges"
    print(io, st)
end

# pretty printing, long version
function Base.show(io::IO, ::MIME"text/plain", ag::AtomGraph)
    st = "AtomGraph $(ag.id) with $(nv(ag.graph)) nodes, $(ne(ag.graph)) edges\n\tatoms: $(ag.elements)"
    print(io, st)
end


"""
    normalized_laplacian(graph)

Compute the normalized graph Laplacian matrix of the input graph, defined as

``I - D^{-1/2} A D^{-1/2}``

where ``A`` is the adjacency matrix and ``D`` is the degree matrix.
"""
function normalized_laplacian(g::G) where {G<:LightGraphs.AbstractGraph}
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

elements(ag::AtomGraph) = ag.elements

# now visualization stuff...

"Get a list of colors to use for graph visualization."
function graph_colors(atno_list, seed_color = colorant"cyan4")
    atom_types = unique(atno_list)
    atom_type_inds = Dict(atom_types[i] => i for i = 1:length(atom_types))
    color_inds = [atom_type_inds[i] for i in atno_list]
    colors = distinguishable_colors(length(atom_types), seed_color)
    return colors[color_inds]
end


"Helper function for sorting because edge ordering isn't preserved when converting to SimpleGraph."
function lt_edge(
    e1::SimpleWeightedGraphs.SimpleWeightedEdge{<:Integer,<:Real},
    e2::SimpleWeightedGraphs.SimpleWeightedEdge{<:Integer,<:Real},
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
    sg = SimpleGraph(adjacency_matrix(ag.graph))
    plt = gplot(
        sg,
        nodefillc = graph_colors(ag.elements),
        nodelabel = ag.elements,
        edgelinewidth = graph_edgewidths(ag),
    )
    display(plt)
end
