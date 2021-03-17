using Test
using LightGraphs
using Serialization
using SimpleWeightedGraphs
using ChemistryFeaturization

@testset "AtomGraph" begin
    # build a silly little triangle graph
    adj = Float32.([0 1 1; 1 0 1; 1 1 0])
    g = SimpleWeightedGraph{Int32}(adj)

    # add an element list that doesn't make sense
    @test_throws AssertionError AtomGraph(g, ["C"])

    # okay, now do it right, start with no features, checking both constructors
    ag = AtomGraph(g, ["C", "C", "C"])
    ag2 = AtomGraph(adj, ["C", "C", "C"])
    @test adjacency_matrix(ag)==adjacency_matrix(ag2)

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
    add_features!(ag3, feature_names, nbins=Int32.([4,3]))

    @test ag.features==ag2.features==ag3.features==Float32.(hcat([vecs["C"] for i in 1:3]...))

    # test that we can add features to multiple graphs
    ag = AtomGraph(g, ["C", "C", "C"])
    ag2 = deepcopy(ag)
    ags = [ag ag2]
    add_features_batch!(ags, featurization)
    
    # and test this one for the other signatures too
    add_features_batch!(ags, vecs, featurization)
    @test ags[1].features[:,1]==Float32.([0;1;0;0;0;1;0])==ags[2].features[:,1]
    # sneakily test that block will always be length 4, i.e. second entry ignored
    add_features_batch!(ags, Symbol.(["X", "Block"]); nbins=[3,3])
    @test ags[1].features[:,2]==Float32.([0;1;0;0;1;0;0])==ags[2].features[:,2]
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
    @test_throws ArgumentError build_graph(joinpath(@__DIR__, "test_data", "nanlaplstruc.cif"))
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
    feature_names = [Symbol("Atomic mass"), :Block]
    featurization = build_featurization(feature_names)
    input_folder = joinpath(@__DIR__, "test_data")
    output_folder=joinpath(@__DIR__, "test_data", "graphs")

    # sneakily run the case wherein output_folder doesn't exist, and overwrite is set to false
    gs = build_graphs_batch(input_folder, featurization, output_folder = output_folder)
    gs2 = build_graphs_batch(input_folder, feature_names, overwrite = true)
    @test repr.(gs)==repr.(gs2) # sneakily testing pretty printing also...

    # test reading from individual files
    g1 = deserialize(joinpath(@__DIR__, "test_data","graphs","mp-195.jls"))
    @test size(g1)==(4,4)
    @test size(g1.features)==(14,4)
    @test g1.id=="mp-195"
    g2 = deserialize(joinpath(@__DIR__, "test_data","graphs","mp-224.jls"))
    @test size(g2)==(6,6)
    @test size(g2.features)==(14,6)
    w = weights(g2)
    @test w[1,2]==w[1,3]==w[1,5]==w[2,4]==w[2,6]==w[3,4]==w[5,6]==0.0
    @test w[3,3]==w[4,4]==w[5,5]==w[6,6]==1.0
    @test g2.elements==["W","W","S","S","S","S"]

    # test read_graphs_batch and alternate syntax of build_graphs_batch
    gs3 = read_graphs_batch(output_folder)
    @test repr(gs2[1])==repr(gs3[1])
    @test repr(gs3[2])==repr(gs[2])

    # test overwrite=true and overwrite=false options
    mtime_initial = map(file -> stat(file).mtime, filter((file) -> isfile(file), readdir(output_folder, join = true)))
    gs4 = build_graphs_batch(input_folder, featurization, output_folder=output_folder, overwrite=true)
    mtime_final = map(file -> stat(file).mtime, filter((file) -> isfile(file), readdir(output_folder, join = true)))
    @test mtime_initial!=mtime_final

    mtime_initial = map(file -> stat(file).mtime, filter((file) -> isfile(file), readdir(output_folder, join = true)))
    gs4 = build_graphs_batch(input_folder, dist_decay_func=exp_decay, output_folder=output_folder, overwrite=false)
    mtime_final = map(file -> stat(file).mtime, filter((file) -> isfile(file), readdir(output_folder, join = true)))
    @test mtime_initial==mtime_final

    rm(output_folder; recursive=true)
end
