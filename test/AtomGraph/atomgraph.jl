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


@testset "save/load" begin
    g = SimpleWeightedGraph{Int32}(Float32.([0 1 1; 1 0 1; 1 1 0]))
    @test_throws AssertionError AtomGraph(g, ["C"])
    fmat = Float32.([1 2 3; 4 5 6; 0 1 0; 9 8 7])
    featurization = [AtomFeat(:feat, true, 2, false, ['a','b']), AtomFeat(:feat2, false, 2, false, [-1,0,1])]
    ag = AtomGraph(g, ["C", "C", "C"], fmat, featurization)
    serialize(abspath(@__DIR__, "../test_data", "testgraph.jls"), ag)
    ag2 = deserialize(abspath(@__DIR__, "../test_data", "testgraph.jls"))
    @test ag2.elements==ag.elements
    @test ag.lapl==ag2.lapl
    @test ag.features==ag2.features
    for i in 1:2
        for field in [:name, :categorical, :num_bins, :logspaced, :vals]
            @test getfield(ag.featurization[i], field)==getfield(ag2.featurization[i], field)
        end
    end
end
