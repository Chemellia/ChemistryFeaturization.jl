using GraphPlot, Colors

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
