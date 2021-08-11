using Test
using ChemistryFeaturization.Utils.GraphBuilding
using Xtals
using Zygote, FiniteDifferences

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

@testset "Graph Building AD tests" begin
  i, j = collect(1:10), collect(1:10)
  dists = Float64.(collect(1:10))

  fd = grad(forward_fdm(2,1),
            (i,j,dist) -> sum(GraphBuilding.weights_cutoff(i,j,dist)),
            i, j, dists)

  gs = gradient(i, j, dists) do i, j, dist
      sum(GraphBuilding.weights_cutoff(i, j, dist))
  end

  @test gs[1] == fill(nothing, 10)
  @test gs[2] == fill(nothing, 10)
  @test gs[3] ≈ fd[3]
end
