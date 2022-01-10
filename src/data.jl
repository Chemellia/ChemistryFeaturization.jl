module Data

using CSV
using DataFrames
using JSON3

export element_data_df, elementfeature_info, valenceshell_conf_df

atom_data_path = joinpath(@__DIR__, "..", "data", "pymatgen_atom_data.csv")
element_data_df = DataFrame(CSV.File(atom_data_path))

elementfeature_info_path = joinpath(@__DIR__, "..", "data", "elementfeature_info.json")
elementfeature_info = JSON3.read(read(elementfeature_info_path, String))

end
