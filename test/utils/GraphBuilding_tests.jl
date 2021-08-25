using Test
using ChemistryFeaturization.Utils.GraphBuilding
using Xtals

@testset "GraphBuilding" begin
    path1 = abspath(@__DIR__, "..", "test_data", "strucs", "mp-195.cif")
    adj, els = build_graph(path1; use_voronoi = true)
    wm_true = [0.0 1.0 1.0 1.0; 1.0 0.0 1.0 1.0; 1.0 1.0 0.0 1.0; 1.0 1.0 1.0 0.0]
    els_true = ["Ho", "Pt", "Pt", "Pt"]

    @test adj == wm_true
    @test els == els_true

    adj, els = build_graph(path1; use_voronoi = false)
    @test adj == wm_true
    @test els == els_true

    # test that we get the same results building from a Crsytal object
    adjc, elsc = build_graph(Crystal(path1))
    @test adjc == wm_true
    @test elsc == els_true
end
