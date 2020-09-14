using Test
using LightGraphs
include("../src/pmg_graphs.jl")
include("../src/atomgraph.jl")

@testset "graph-building" begin
    wm, atoms = build_graph(joinpath(@__DIR__, "./test_data/mp-195.cif"))
    wm_true = [0.0 1.0 1.0 1.0; 1.0 0.0 1.0 1.0; 1.0 1.0 0.0 1.0; 1.0 1.0 1.0 0.0]
    @test wm == wm_true
    @test atoms == ["Ho", "Pt", "Pt", "Pt"]
    wm, atoms = build_graph(joinpath(@__DIR__, "./test_data/mp-195.cif"); use_voronoi=false)
    @test wm == wm_true
    @test atoms == ["Ho", "Pt", "Pt", "Pt"]
end

@testset "AtomGraph" begin
    # build a silly little graph and make sure fcns evaluate as expected
    g = SimpleWeightedGraph(Float32.([0 1 1; 1 0 1; 1 1 0]))

    # add an element list that doesn't make sense

    # okay, now do it right, stat with no features
    ag = AtomGraph(g, ["C", "C", "C"])

    # check stuff


    # add some bad features


    # okay, doing it right again...

end