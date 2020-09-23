# storing this here for now
using JSON3, ChemistryFeaturization, SimpleWeightedGraphs, StructTypes
g = SimpleWeightedGraph{Int32}(Float32.([0 1 1; 1 0 1; 1 1 0]))
fmat = Float32.([1 2 3; 4 5 6])
featurization = [AtomFeat(:feat, true, 2, false, ['a','b'])]
ag = AtomGraph(g, ["C", "C", "C"], fmat, featurization)

StructTypes.StructType(::Type{AtomGraph})=StructTypes.Struct()
StructTypes.StructType(::Type{AtomFeat})=StructTypes.Struct()
StructTypes.StructType(::Type{<:SimpleWeightedGraph})=StructTypes.Struct()
#StructTypes.StructType(::Type{<:SparseArrays.SparseMatrixCSC})=StructTypes.ArrayType()
StructTypes.StructType(::Type{<:SparseArrays.SparseMatrixCSC})=StructTypes.DictType()

s = JSON3.write(ag)
ag2 = JSON3.read(s, AtomGraph)