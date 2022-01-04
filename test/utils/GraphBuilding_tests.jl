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

  function test_fd(i, j, dist)
      fd = grad(forward_fdm(2,1),
                (i,j,dist) -> sum(GraphBuilding.weights_cutoff(i,j,dist)),
                i, j, dist)

      gs = gradient(i, j, dist) do i, j, dist
          sum(GraphBuilding.weights_cutoff(i, j, dist))
      end

      @test gs[1] == nothing # fill(nothing, length(i))
      @test gs[2] == nothing # fill(nothing, length(j))
      t = isapprox(gs[3], fd[3])
      @test t
  end

  # test with non-overlapping indices
  test_fd(collect(1:10), collect(1:10), Float64.(collect(1:10)))
  # test with overlapping indices
  test_fd(rand(1:10, 100), rand(1:10, 100), rand(100))
end
