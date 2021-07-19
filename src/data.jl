module Data

using CSV
using DataFrames
using JSON

export atom_data_df, feature_info

atom_data_path = joinpath(@__DIR__, "..", "data", "pymatgen_atom_data.csv")
atom_data_df = DataFrame(CSV.File(atom_data_path))

feature_info_path = joinpath(@__DIR__, "..", "data", "feature_info.json")
feature_info = JSON.parsefile(feature_info_path)

end
