module OrbitalFeatureUtils

using DataFrames

# Helper functions used for encoding
#=
return the Index and Value arrays (can be used for building a sparse matrix on their own)
corresponding to a single element given its valence shell configuration (regex-extracted usually)
=#
function _orbitalsparse(valence_shell_config)
    I = map(e -> _orbitalindex(e[1:2]), valence_shell_config)
    V = map(e -> parse(Int, e[3]), valence_shell_config)
    return I, V
end

# faster and more memory efficient than using a dictionary
function _orbitalindex(orbital::SubString{String})
    if orbital == "1s"
        return 1
    elseif orbital == "2s"
        return 2
    elseif orbital == "2p"
        return 3
    elseif orbital == "3s"
        return 4
    elseif orbital == "3p"
        return 5
    elseif orbital == "4s"
        return 6
    elseif orbital == "3d"
        return 7
    elseif orbital == "4p"
        return 8
    elseif orbital == "5s"
        return 9
    elseif orbital == "4d"
        return 10
    elseif orbital == "5p"
        return 11
    elseif orbital == "6s"
        return 12
    elseif orbital == "4f"
        return 13
    elseif orbital == "5d"
        return 14
    elseif orbital == "6p"
        return 15
    elseif orbital == "7s"
        return 16
    elseif orbital == "5f"
        return 17
    elseif orbital == "6d"
        return 18
    elseif orbital == "7p"
        return 19
    elseif orbital == "8s"
        return 20
    elseif orbital == "5g"
        return 21
    elseif orbital == "6f"
        return 22
    elseif orbital == "7d"
        return 23
    elseif orbital == "8p"
        return 24
    else # orbital == "9s"
        return 25
    end
end


# using regex, get the valence shell configuration
function _orbitalregex(df::DataFrame, element)
    config_df_row = filter(:Symbol => x -> x == element, df; view = true)[
        !,
        "Valence Shell Configuration",
    ][1]
    orbital_regex = r"[1-9][s|p|d|f|g][1-9]" # regex which will match individual orbitals
    orbital_config = map(o -> o.match, collect(eachmatch(orbital_regex, config_df_row)))
    return orbital_config
end

# should this be moved into something like utils/orbitalfeatureutils.jl?
function valence_shell_config(df::DataFrame, element)
    valence_shell_config = _orbitalregex(df, element)
    return _orbitalsparse(valence_shell_config)
end


# Helper functions used for decoding

function _indexorbital(orbital::Int16)
    if orbital == 1
        return "1s"
    elseif orbital == 2
        return "2s"
    elseif orbital == 3
        return "2p"
    elseif orbital == 4
        return "3s"
    elseif orbital == 5
        return "3p"
    elseif orbital == 6
        return "4s"
    elseif orbital == 7
        return "3d"
    elseif orbital == 8
        return "4p"
    elseif orbital == 9
        return "5s"
    elseif orbital == 10
        return "4d"
    elseif orbital == 11
        return "5p"
    elseif orbital == 12
        return "6s"
    elseif orbital == 13
        return "4f"
    elseif orbital == 14
        return "5d"
    elseif orbital == 15
        return "6p"
    elseif orbital == 16
        return "7s"
    elseif orbital == 17
        return "5f"
    elseif orbital == 18
        return "6d"
    elseif orbital == 19
        return "7p"
    elseif orbital == 20
        return "8s"
    elseif orbital == 21
        return "5g"
    elseif orbital == 22
        return "6f"
    elseif orbital == 23
        return "7d"
    elseif orbital == 24
        return "8p"
    else
        return "9s"
    end
end

end
