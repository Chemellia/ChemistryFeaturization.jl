using Test
using LightGraphs
using Serialization
using SimpleWeightedGraphs
using ChemistryFeaturization

feature_names = [Symbol("Atomic mass"), :Block]
featurization = build_featurization(feature_names)
input_folder = abspath(@__DIR__, "../test_data")
output_folder=abspath(@__DIR__, "../test_data", "graphs")


@testset "initial batch processing" begin
    # sneakily run the case wherein output_folder doesn't exist, and overwrite is set to false
    gs = build_graphs_batch(input_folder, featurization, output_folder = output_folder)
    gs2 = build_graphs_batch(input_folder, feature_names, overwrite = true)
    @test repr.(gs)==repr.(gs2) # sneakily testing pretty printing also...

    rm(output_folder; recursive=true)

    gs = build_graphs_batch(input_folder, featurization, output_folder = output_folder)
    gs2 = build_graphs_batch(input_folder, feature_names, output_folder = output_folder, overwrite = false)
    @test repr.(gs)==repr.(gs2) # sneakily testing pretty printing also...
end


@testset "deserialization" begin
    # test reading from individual files
    g1 = deserialize(abspath(@__DIR__, "../test_data","graphs","mp-195.jls"))
    @test size(g1)==(4,4)
    @test size(g1.features)==(14,4)
    @test g1.id=="mp-195"
    g2 = deserialize(abspath(@__DIR__, "../test_data","graphs","mp-224.jls"))
    @test size(g2)==(6,6)
    @test size(g2.features)==(14,6)
    w = weights(g2)
    @test w[1,2]==w[1,3]==w[1,5]==w[2,4]==w[2,6]==w[3,4]==w[5,6]==0.0
    @test w[3,3]==w[4,4]==w[5,5]==w[6,6]==1.0
    @test g2.elements==["W","W","S","S","S","S"]

end


@testset "read graphs" begin
    # test read_graphs_batch and alternate syntax of build_graphs_batch
    gs = build_graphs_batch(input_folder, featurization, output_folder = output_folder)
    gs2 = build_graphs_batch(input_folder, feature_names, output_folder = output_folder, overwrite = false)
    gs3 = read_graphs_batch(output_folder)
    @test repr(gs2[1])==repr(gs3[1])
    @test repr(gs3[2])==repr(gs[2])
end


@testset "batch processing - overwrite" begin
    # test overwrite=true option
    mtime_initial = map(file -> stat(file).mtime, filter((file) -> isfile(file), readdir(output_folder, join = true)))
    gs4 = build_graphs_batch(input_folder, featurization, output_folder=output_folder, overwrite=true)
    mtime_final = map(file -> stat(file).mtime, filter((file) -> isfile(file), readdir(output_folder, join = true)))
    @test mtime_initial!=mtime_final

    # test overwrite=false option 
    mtime_initial = map(file -> stat(file).mtime, filter((file) -> isfile(file), readdir(output_folder, join = true)))
    gs4 = build_graphs_batch(input_folder, dist_decay_func=exp_decay, output_folder=output_folder, overwrite=false)
    mtime_final = map(file -> stat(file).mtime, filter((file) -> isfile(file), readdir(output_folder, join = true)))
    @test mtime_initial==mtime_final
    rm(output_folder; recursive=true)
end
