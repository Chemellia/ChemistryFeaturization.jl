using Test
using LightGraphs
using SimpleWeightedGraphs
using Serialization
using ..ChemistryFeaturization.Atoms

# NB: featurizing graphs is tested in ElementFeatureDescriptor and GraphNodeFeaturization tests

@testset "AtomGraph" begin

    wm_mp195 = [0.0 1.0 1.0 1.0; 1.0 0.0 1.0 1.0; 1.0 1.0 0.0 1.0; 1.0 1.0 1.0 0.0]
    lapl_mp195 = [
        1.0 -0.3333333 -0.3333333 -0.3333333
        -0.3333333 1.0 -0.3333333 -0.3333333
        -0.3333333 -0.3333333 1.0 -0.3333333
        -0.3333333 -0.3333333 -0.3333333 1.0
    ]
    els_mp195 = ["Ho", "Pt", "Pt", "Pt"]

    @testset "construct object" begin
        # build a cute little triangle graph
        adj = Float32.([0 1 1; 1 0 1; 1 1 0])
        g = SimpleWeightedGraph(adj)

        # add an element list that doesn't make sense
        @test_throws AssertionError AtomGraph(g, ["C"])

        # okay, now do it right, start with no features, checking both constructors
        ag = AtomGraph(g, ["C", "C", "C"])
        ag2 = AtomGraph(adj, ["C", "C", "C"])
        @test adjacency_matrix(ag.graph) == adjacency_matrix(ag2.graph)

        # and now from a file...
        ag = AtomGraph(abspath(@__DIR__, "..", "test_data", "strucs", "mp-195.cif"))

        @test weights(ag.graph) == wm_mp195
        @test all(isapprox.(ag.laplacian, lapl_mp195, atol = 1e-7))
        @test elements(ag) == els_mp195

        # test that warning is thrown for NaNs in laplacian
        @test_throws ArgumentError AtomGraph(
            build_graph(
                abspath(@__DIR__, "..", "test_data", "strucs", "nanlaplstruc.cif"),
            )...,
        )
    end

    @testset "save/load" begin
        g = SimpleWeightedGraph(Float32.([0 1 1; 1 0 1; 1 1 0]))
        ag = AtomGraph(g, ["C", "C", "C"])
        serialize(abspath(@__DIR__, "..", "test_data", "strucs", "testgraph.jls"), ag)
        ag2 = deserialize(abspath(@__DIR__, "..", "test_data", "strucs", "testgraph.jls"))
        @test adjacency_matrix(ag.graph) == adjacency_matrix(ag2.graph)
        @test elements(ag) == ag2.elements
        @test Atoms.normalized_laplacian(ag) == ag2.laplacian
    end

    @testset "batch processing" begin
        file_list = readdir(abspath(@__DIR__, "..", "test_data", "strucs"), join = true)
        ags = AtomGraph.(file_list)
        @test length(collect(skipmissing(ags))) == length(ags) - 1

        ids = map(f -> splitpath(f)[end], file_list)
        ags2 = collect(skipmissing(AtomGraph.(file_list, ids)))
        mp195s = [ag for ag in skipmissing(ags2) if contains(ag.id, "mp-195")]

        for ag in mp195s
            @test elements(ag) == els_mp195
            @test weights(ag.graph) == wm_mp195
            @test all(isapprox.(ag.laplacian, lapl_mp195, atol = 1e-7))
        end

    end

    @testset "visualize" begin
        # test lt_edge
        adj = Float32.([0 1 2; 1 0 1; 2 1 0])
        g = SimpleWeightedGraph{Int32}(adj)
        e_adj = [e for e in edges(g)]

        @test ChemistryFeaturization.Atoms.lt_edge(e_adj[2], e_adj[3]) == true # e1.src < e2.src
        @test ChemistryFeaturization.Atoms.lt_edge(e_adj[1], e_adj[2]) == true # e1.dest < e2.dest
        @test ChemistryFeaturization.Atoms.lt_edge(e_adj[3], e_adj[2]) == false
        adj = Float32.([0 1 1; 1 0 1; 1 1 0])
        g = SimpleWeightedGraph{Int32}(adj)
        ag = AtomGraph(g, ["C", "C", "C"])
        @test ChemistryFeaturization.Atoms.graph_edgewidths(ag) == [1.0, 1.0, 1.0]

        # test graph_colors
        ag = AtomGraph(g, ["O", "C", "O"])
        c1, c2, c3 = ChemistryFeaturization.Atoms.graph_colors(elements(ag))
        @test c1 == c3 != c2

        # TODO: figure out a way to test visualize_graph itself?
    end

end
