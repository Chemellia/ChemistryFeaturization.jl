module Data

using CSV
using DataFrames
using JSON3

export atom_data_df, feature_info, valenceshell_conf_df

atom_data_path = joinpath(@__DIR__, "..", "data", "pymatgen_atom_data.csv")
atom_data_df = DataFrame(CSV.File(atom_data_path))

feature_info_path = joinpath(@__DIR__, "..", "data", "feature_info.json")
feature_info = JSON3.read(read(feature_info_path, String))

valenceshell_conf_df = select(
    atom_data_df,
    :Symbol,
    "Electronic Structure" =>
        (
            x -> replace.(x, r"\[(.+)\]\.(?<valence>\w+)" => s"\g<valence>")
        ) => "Electronic Structure",
)

end
