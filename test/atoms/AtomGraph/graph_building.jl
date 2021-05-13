using Test
using LightGraphs
using Serialization
using SimpleWeightedGraphs
using ChemistryFeaturization


@testset "graph-building" begin
    ag = build_graph(abspath(@__DIR__, "../test_data", "mp-195.cif"), use_voronoi=true)
    wm_true = [0.0 1.0 1.0 1.0; 1.0 0.0 1.0 1.0; 1.0 1.0 0.0 1.0; 1.0 1.0 1.0 0.0]
    @test weights(ag) == wm_true
    @test ag.elements == ["Ho", "Pt", "Pt", "Pt"]
    ag = build_graph(abspath(@__DIR__, "../test_data", "mp-195.cif"); use_voronoi=false)
    @test weights(ag) == wm_true
    @test ag.elements == ["Ho", "Pt", "Pt", "Pt"]

    # tests for some other file formats
    graphs = AtomGraph[]
    for fp in ["mp-195.poscar", "mp-195.traj", "mp-195.xyz"]
        push!(graphs, build_graph(abspath(@__DIR__, "../test_data", fp)))
    end
    for g in graphs
        @test ag.elements == g.elements
        @test weights(ag) == weights(g)
    end
end


@testset "non-periodic system graph-building" begin
    # test for nonperiodic system
    @test_logs (:warn, "Voronoi edge weights are not supported if any direction in the structure is nonperiodic. Using cutoff weights method...") build_graph(abspath(@__DIR__,"../test_data", "methane.xyz"), use_voronoi=true)
    methane = build_graph(abspath(@__DIR__,"../test_data", "methane.xyz"))
    @test all(isapprox.(weights(methane)[2:5,1], 1.0, atol=1e-4))
    @test all(isapprox.(weights(methane)[3:2,2], 0.375, atol=1e-5))
end


@testset "graph-building - NaNs in Laplacian" begin
    # test that warning is thrown for NaNs in laplacian
    @test_throws ArgumentError build_graph(abspath(@__DIR__, "../test_data", "nanlaplstruc.cif"))
end
