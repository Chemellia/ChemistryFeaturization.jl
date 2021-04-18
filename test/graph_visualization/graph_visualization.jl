using Test
using LightGraphs
using Serialization
using SimpleWeightedGraphs
using ChemistryFeaturization

# test lt_edge
adj = Float32.([0 1 2; 1 0 1; 2 1 0])
g = SimpleWeightedGraph{Int32}(adj)
e_adj = [e for e in edges(g)]

@test ChemistryFeaturization.lt_edge(e_adj[2], e_adj[3])==true # e1.src < e2.src
@test ChemistryFeaturization.lt_edge(e_adj[1], e_adj[2])==true # e1.dest < e2.dest
@test ChemistryFeaturization.lt_edge(e_adj[3], e_adj[2])==false

# test graph_edgewidths
adj = Float32.([0 1 1; 1 0 1; 1 1 0])
g = SimpleWeightedGraph{Int32}(adj)
ag = AtomGraph(g, ["C", "C", "C"])
@test ChemistryFeaturization.graph_edgewidths(ag)==[1.0, 1.0, 1.0]

# test graph_colors
ag = AtomGraph(g, ["O", "C", "O"])
c1, c2, c3 = ChemistryFeaturization.graph_colors(ag.elements)
@test c1 == c3 != c2

# visualize_graph yet to be tested
