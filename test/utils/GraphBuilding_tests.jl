using Test
using ChemistryFeaturization.Utils.GraphBuilding
using Zygote, FiniteDifferences

@testset "GraphBuilding" begin
    adj, els = build_graph(
        abspath(@__DIR__, "..", "test_data", "strucs", "mp-195.cif"),
        use_voronoi = true,
    )
    wm_true = [0.0 1.0 1.0 1.0; 1.0 0.0 1.0 1.0; 1.0 1.0 0.0 1.0; 1.0 1.0 1.0 0.0]
    @test adj == wm_true
    @test els == ["Ho", "Pt", "Pt", "Pt"]

    adj, els = build_graph(
        abspath(@__DIR__, "..", "test_data", "strucs", "mp-195.cif");
        use_voronoi = false,
    )
    @test adj == wm_true
    @test els == ["Ho", "Pt", "Pt", "Pt"]

    # tests for some other file formats
    info = Tuple{Matrix,Vector{String}}[]
    for fp in ["mp-195.poscar", "mp-195.traj", "mp-195.xyz"]
        push!(info, build_graph(abspath(@__DIR__, "..", "test_data", "strucs", fp)))
    end
    for t in info
        @test t[1] == wm_true
        @test t[2] == els
    end

    # test for nonperiodic system
    @test_logs (
        :warn,
        "Voronoi edge weights are not supported if any direction in the structure is nonperiodic. Using cutoff weights method...",
    ) build_graph(
        abspath(@__DIR__, "..", "test_data", "strucs", "methane.xyz"),
        use_voronoi = true,
    )
    adj, els = build_graph(abspath(@__DIR__, "..", "test_data", "strucs", "methane.xyz"))
    @test all(isapprox.(adj[2:5, 1], 1.0, atol = 1e-4))
    @test all(isapprox.(adj[3:2, 2], 0.375, atol = 1e-5))
end

@testset "Graph Building AD tests" begin

  function test_fd(i, j, dist)
      fd = grad(forward_fdm(2,1),
                (i,j,dist) -> sum(GraphBuilding.weights_cutoff(i,j,dist)),
                i, j, dist)

      gs = gradient(i, j, dists) do i, j, dist
          sum(GraphBuilding.weights_cutoff(i, j, dist))
      end

      @test gs[1] == fill(nothing, length(i))
      @test gs[2] == fill(nothing, length(j))
      @test gs[3] ≈ fd[3]
  end

  # test with non-overlapping indices
  test_fd(collect(1:10), collect(1:10), Float64.(collect(1:10)))
  # test with overlapping indices
  test_fd(rand(1:10, 100), rand(1:10, 100), rand(100))
end
