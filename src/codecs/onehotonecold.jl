using Flux: onecold

include("onehotonecold_utils.jl")

"""
    OneHotOneCold(categorical, bins)

AbstractCodec type which uses a dummy variable (as defined in statistical literature), i.e., which employs
one-hot encoding and a one-cold decoding scheme.
"""
struct OneHotOneCold <: AbstractCodec
    categorical::Bool
    bins::Vector
    function OneHotOneCold(categorical, bins)
        @assert categorical || eltype(bins) <: Number "Your bins aren't numbers...are you sure you didn't mean for this codec to be categorical?"
        new(categorical, bins)
    end
end

output_shape(ohoc::OneHotOneCold) =
    ohoc.categorical ? length(ohoc.bins) : length(ohoc.bins) - 1

"A flexible version of Flux.onehot that can handle both categorical and continuous-valued encoding."
function encode(val, ohoc::OneHotOneCold)
    local bin_index, onehot_vec
    if ohoc.categorical
        onehot_vec = [0.0 for i = 1:length(ohoc.bins)]
        bin_index = findfirst(isequal(val), ohoc.bins)
    else
        @assert ohoc.bins[1] <= val <= ohoc.bins[end] "The value $val is outside the range of bins $ohoc.bins"
        onehot_vec = [0.0 for i = 1:(length(ohoc.bins)-1)]
        bin_index = searchsorted(ohoc.bins, val).stop
        if bin_index == length(ohoc.bins) # got the max value
            bin_index = bin_index - 1
        elseif isapprox(val, ohoc.bins[1]) # sometimes will get 0 if this doesn't get checked
            bin_index = 1
        end
    end
    onehot_vec[bin_index] = 1.0
    return onehot_vec
end

encode(::Missing, ohoc::OneHotOneCold) = Vector{Missing}(missing, output_shape(ohoc))

encode(vals::Vector, ohoc::OneHotOneCold) = encode.(vals, Ref(ohoc))

function encode(vals::Array, ohoc::OneHotOneCold)
    # output is one dimension larger
    output = Array{Union{Bool,Missing}}(missing, size(vals)..., output_shape(ohoc))
    for ind in CartesianIndices(vals)
        output[ind, :] = encode(vals[ind], ohoc)
    end
    return output
end

function decode(encoded::Vector, ohoc::OneHotOneCold)
    @assert length(encoded) == output_shape(ohoc)
    local decoded
    if ohoc.categorical # return value
        decoded = onecold(encoded, ohoc.bins)
    else # return bounds (TODO: should this be a tuple or a vector..? I like tuple for distinguishing from encoded vectors, but it doesn't play so nice with broadcasting...)
        decoded = (onecold(encoded, ohoc.bins[1:end-1]), onecold(encoded, ohoc.bins[2:end]))
    end
    return decoded
end

# this is separate from the Array case because at the moment there are different conventions of index ordering for ndims == 2 vs. ndims > 2...oops, we should probably change that, but it will require changes in AtomicGraphNets as well
function decode(encoded::Matrix, ohoc::OneHotOneCold)
    @assert size(encoded)[1] == output_shape(ohoc)
    local decoded
    decoded_length = size(encoded)[2]
    decoded_eltype = eltype(ohoc.bins)
    if ohoc.categorical
        decoded = Vector{Union{decoded_eltype,Missing}}(missing, decoded_length)
    else
        decoded = Vector{Union{Tuple{decoded_eltype,decoded_eltype},Missing}}(
            missing,
            decoded_length,
        )
    end
    for i = 1:decoded_length
        vec = encoded[:, i] # this is the key difference from the Array dispatch below
        if !all(ismissing.(vec))
            decoded[i] = decode(vec, ohoc)
        end
    end
    return decoded
end

function decode(encoded::Array, ohoc::OneHotOneCold)
    @assert size(encoded)[end] == output_shape(ohoc)
    local decoded
    decoded_size = size(encoded)[1:end-1]
    decoded_eltype = eltype(ohoc.bins)
    if ohoc.categorical
        decoded = Array{Union{decoded_eltype,Missing}}(missing, decoded_size...)
    else
        decoded = Array{Union{Tuple{decoded_eltype,decoded_eltype},Missing}}(
            missing,
            decoded_size...,
        )
    end
    for ind in CartesianIndices(decoded)
        vec = encoded[ind, :]
        if !all(ismissing.(vec))
            decoded[ind] = decode(vec, ohoc)
        end
    end
    return decoded
end
