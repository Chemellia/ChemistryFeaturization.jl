using Test
using LightGraphs
using Serialization
using SimpleWeightedGraphs
using ChemistryFeaturization

# featurizing graphs is tested in feature object and featurization tests

@testset "AtomGraph" begin

    @testset "construct_object" begin
        # build a silly little triangle graph
        adj = Float32.([0 1 1; 1 0 1; 1 1 0])
        g = SimpleWeightedGraph{Int32}(adj)

        # add an element list that doesn't make sense
        @test_throws AssertionError AtomGraph(g, ["C"])

        # okay, now do it right, start with no features, checking both constructors
        ag = AtomGraph(g, ["C", "C", "C"])
        ag2 = AtomGraph(adj, ["C", "C", "C"])
        @test adjacency_matrix(ag.graph)==adjacency_matrix(ag2.graph)
    end

    @testset "save/load" begin
        g = SimpleWeightedGraph{Int32}(Float32.([0 1 1; 1 0 1; 1 1 0]))
        ag = AtomGraph(g, ["C", "C", "C"])
        serialize(abspath(@__DIR__, "../test_data", "testgraph.jls"), ag)
        ag2 = deserialize(abspath(@__DIR__, "../test_data", "testgraph.jls"))
        @test adjacency_matrix(ag.graph).==adjacency_matrix(ag2.graph)
        @test ag2.elements==ag.elements
        @test ag.laplacian==ag2.laplacian
    end

    @testset "visualize" begin
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
    end

    @testset "batch_processing" begin

    end

    @testset "build_adjacency_matrix" begin

    end

end

