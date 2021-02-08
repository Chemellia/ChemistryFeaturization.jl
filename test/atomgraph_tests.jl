using Test
using LightGraphs
using Serialization
include("../src/pmg_graphs.jl")
include("../src/atomgraph.jl")

@testset "AtomGraph" begin
    # build a silly little triangle graph
    g = SimpleWeightedGraph{Int32}(Float32.([0 1 1; 1 0 1; 1 1 0]))

    # add an element list that doesn't make sense
    @test_throws AssertionError AtomGraph(g, ["C"])

    # okay, now do it right, start with no features
    ag = AtomGraph(g, ["C", "C", "C"])

    # check LightGraphs fcns
    @test eltype(ag)==Int32
    @test edgetype(ag)==SimpleWeightedEdge{Int32,Float32}
    @test ne(ag)==3
    @test nv(ag)==3
    @test !is_directed(ag)
    # not sure the best way to test the ones that return iterators, e.g. edges
    @test outneighbors(ag,1)==inneighbors(ag,1)==[2,3]
    @test has_vertex(ag,1)
    @test !has_vertex(ag,4)
    @test has_edge(ag,1,2)

    # add some features
    bad_fmat = Float32.([1 2; 3 4])
    good_fmat = Float32.([1 2 3; 4 5 6])
    featurization = [AtomFeat(:feat, true, 2, false, ['a','b'])]
    @test_throws AssertionError add_features!(ag, bad_fmat, featurization)
    add_features!(ag, good_fmat, featurization)
    @test ag.features==good_fmat

    # tests for other signatures of add_features! where feature vectors are built automatically
    ag = AtomGraph(g, ["C", "C", "C"])
    ag2 = deepcopy(ag)
    ag3 = deepcopy(ag)
    feature_names = [:Block, :X]
    vecs, featurization = make_feature_vectors(feature_names, nbins=[4,3])

    add_features!(ag, vecs, featurization)
    add_features!(ag2, featurization)
    add_features!(ag3, feature_names, Int32.([4,3]))

    @test ag.features==ag2.features==ag3.features==Float32.(hcat([vecs["C"] for i in 1:3]...))
end

@testset "graph-building" begin
    ag = build_graph(joinpath(@__DIR__, "test_data", "mp-195.cif"), use_voronoi=true)
    wm_true = [0.0 1.0 1.0 1.0; 1.0 0.0 1.0 1.0; 1.0 1.0 0.0 1.0; 1.0 1.0 1.0 0.0]
    @test weights(ag) == wm_true
    @test ag.elements == ["Ho", "Pt", "Pt", "Pt"]
    ag = build_graph(joinpath(@__DIR__, "test_data", "mp-195.cif"); use_voronoi=false)
    @test weights(ag) == wm_true
    @test ag.elements == ["Ho", "Pt", "Pt", "Pt"]

    # tests for some other file formats
    graphs = AtomGraph[]
    for fp in ["mp-195.poscar", "mp-195.traj", "mp-195.xyz"]
        push!(graphs, build_graph(joinpath(@__DIR__, "test_data", fp)))
    end
    for g in graphs
        @test ag.elements == g.elements
        @test weights(ag) == weights(g)
    end

    # test for nonperiodic system
    @test_logs (:warn, "Voronoi edge weights are not supported if any direction in the structure is nonperiodic. Using cutoff weights method...") build_graph(joinpath(@__DIR__,"test_data", "methane.xyz"), use_voronoi=true)
    methane = build_graph(joinpath(@__DIR__,"test_data", "methane.xyz"))
    @test all(isapprox.(weights(methane)[2:5,1], 1.0, atol=1e-4))
    @test all(isapprox.(weights(methane)[3:2,2], 0.375, atol=1e-5))

    # test that warning is thrown for NaNs in laplacian
    @test_throws AssertionError build_graph(joinpath(@__DIR__, "test_data", "nanlaplstruc.cif"))

end

@testset "save/load" begin
    g = SimpleWeightedGraph{Int32}(Float32.([0 1 1; 1 0 1; 1 1 0]))
    @test_throws AssertionError AtomGraph(g, ["C"])
    fmat = Float32.([1 2 3; 4 5 6; 0 1 0; 9 8 7])
    featurization = [AtomFeat(:feat, true, 2, false, ['a','b']), AtomFeat(:feat2, false, 2, false, [-1,0,1])]
    ag = AtomGraph(g, ["C", "C", "C"], fmat, featurization)
    serialize(joinpath(@__DIR__, "test_data", "testgraph.jls"), ag)
    ag2 = deserialize(joinpath(@__DIR__, "test_data", "testgraph.jls"))
    @test ag2.elements==ag.elements
    @test ag.lapl==ag2.lapl
    @test ag.features==ag2.features
    for i in 1:2
        for field in [:name, :categorical, :num_bins, :logspaced, :vals]
            @test getfield(ag.featurization[i], field)==getfield(ag2.featurization[i], field)
        end
    end
end

@testset "batch processing" begin
    featurization = build_atom_feats([Symbol("Atomic mass"), :Block])
    build_graphs_batch(joinpath(@__DIR__, "test_data"), joinpath(@__DIR__, "test_data", "graphs"), featurization)

    # try to figure out what Windows test is doing that it can't find the file
    println(readdir(joinpath(@__DIR__, "test_data", "graphs")

    g1 = deserialize(joinpath(@__DIR__, "test_data","graphs","mp-195.jls"))
    @test size(g1)==(4,4)
    @test size(g1.features)==(14,4)
    g2 = deserialize(joinpath(@__DIR__, "test_data","graphs","mp-224.jls"))
    @test size(g2)==(6,6)
    @test size(g2.features)==(14,6)
    w = weights(g2)
    @test w[1,2]==w[1,3]==w[1,5]==w[2,4]==w[2,6]==w[3,4]==w[5,6]==0.0
    @test w[3,3]==w[4,4]==w[5,5]==w[6,6]==1.0
    @test g2.elements==["W","W","S","S","S","S"]

    # test read_graphs_batch
    gs = read_graphs_batch(joinpath(@__DIR__, "test_data", "graphs"))
    @test length(gs)>=2
end
