using ChemistryFeaturization
using SimpleWeightedGraphs
using JSON3
using JLD2
using BenchmarkTools

# JLD2 version
function saveload_JLD2(ag)
    @save "./test_data/testgraph.jld2" ag
    @load "./test_data/testgraph.jld2" ag
end

"""
# JSON version...haven't gotten this working yet
function saveload_JSON(ag)
    open("./test_data/testgraph.json", "w") do f
        JSON.print(f, ag)
    end
    d = JSON.parsefile("./test_data/testgraph.json"; inttype=Int32)
    ag = AtomGraph(SimpleWeightedGraph{Int32,Float32}(Float32.(d["graph"]["weights"])), String.(d["elements"]), Float32.(d["lapl"]), )
end
"""
# first, a silly tiny example
g = SimpleWeightedGraph{Int32}(Float32.([0 1 1; 1 0 1; 1 1 0]))
fmat = Float32.([1 2 3; 4 5 6])
featurization = [AtomFeat(:feat, true, 2, false, ['a','b'])]
ag = AtomGraph(g, ["C", "C", "C"], fmat, featurization)



