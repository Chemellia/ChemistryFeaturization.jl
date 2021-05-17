using Test

@testset "GraphBuilding" begin
    adj, els = build_graph(abspath(@__DIR__, "../test_data", "mp-195.cif"), use_voronoi=true)
    wm_true = [0.0 1.0 1.0 1.0; 1.0 0.0 1.0 1.0; 1.0 1.0 0.0 1.0; 1.0 1.0 1.0 0.0]
    @test adj == wm_true
    @test els == ["Ho", "Pt", "Pt", "Pt"]

    adj, els = build_graph(abspath(@__DIR__, "../test_data", "mp-195.cif"); use_voronoi=false)
    @test adj == wm_true
    @test els == ["Ho", "Pt", "Pt", "Pt"]

    # tests for some other file formats
    info = Tuple{Matrix,Vector{String}}[]
    for fp in ["mp-195.poscar", "mp-195.traj", "mp-195.xyz"]
        push!(info, build_graph(abspath(@__DIR__, "../test_data", fp)))
    end
    for t in info
        @test t[1] == wm_true
        @test t[2] == els
    end

    # test for nonperiodic system
    @test_logs (:warn, "Voronoi edge weights are not supported if any direction in the structure is nonperiodic. Using cutoff weights method...") build_graph(abspath(@__DIR__,"../test_data", "methane.xyz"), use_voronoi=true)
    adj, els = build_graph(abspath(@__DIR__,"../test_data", "methane.xyz"))
    @test all(isapprox.(adj[2:5,1], 1.0, atol=1e-4))
    @test all(isapprox.(adj[3:2,2], 0.375, atol=1e-5))
end