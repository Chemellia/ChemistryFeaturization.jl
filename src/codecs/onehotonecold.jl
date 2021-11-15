using ..ChemistryFeaturization.AbstractType: AbstractCodec
import ..ChemistryFeaturization.encode
import ..ChemistryFeaturization.decode
using Flux: onecold

export encode, decode, output_shape

"""
    OneHotOneCold(categorical, bins)

AbstractCodec type which uses a dummy variable (as defined in statistical literature), i.e., which employs
one-hot encoding and a one-cold decoding scheme.
"""
struct OneHotOneCold <: AbstractCodec
    categorical::Bool
    bins::Vector
end

output_shape(ohoc::OneHotOneCold) =
    ohoc.categorical ? length(ohoc.bins) : length(ohoc.bins) - 1

"A flexible version of Flux.onehot that can handle both categorical and continuous-valued encoding."
function encode(ohoc::OneHotOneCold, val)
    local bin_index, onehot_vec
    if ohoc.categorical
        onehot_vec = [0.0 for i = 1:length(ohoc.bins)]
        bin_index = findfirst(isequal(val), ohoc.bins)
    else
        @assert eltype(ohoc.bins) <: Number "Your bins aren't numbers...are you sure you didn't mean for this feature to be categorical?"
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

encode(ohoc::OneHotOneCold, ::Missing) = Vector{Missing}(missing, output_shape(ohoc))

encode(ohoc::OneHotOneCold, vals::Vector) = encode.(Ref(ohoc), vals)

function encode(ohoc::OneHotOneCold, vals::Array)
    # output is one dimension larger
    output = Array{Union{Bool,Missing}}(missing, size(vals)..., output_shape(ohoc))
    for ind in CartesianIndices(vals)
        output[ind,:] = encode(ohoc, vals[ind])
    end
    return output
end

function decode(ohoc::OneHotOneCold, encoded::Vector)
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
function decode(ohoc::OneHotOneCold, encoded::Matrix)
    @assert size(encoded)[1] == output_shape(ohoc)
    local decoded
    decoded_length = size(encoded)[2]
    decoded_eltype = eltype(ohoc.bins)
    if ohoc.categorical
        decoded = Vector{Union{decoded_eltype,Missing}}(missing, decoded_length)
    else
        decoded = Array{Union{Tuple{decoded_eltype,decoded_eltype},Missing}}(missing, length)
    end
    for i in 1:decoded_length
        vec = encoded[:,i] # this is the key difference from the Array dispatch below
        if !all(ismissing.(vec))
            decoded[i] = decode(ohoc, vec)
        end
    end
    return decoded
end

function decode(ohoc::OneHotOneCold, encoded::Array)
    @assert size(encoded)[end] == output_shape(ohoc)
    local decoded
    decoded_size = size(encoded)[1:end-1]
    decoded_eltype = eltype(ohoc.bins)
    if ohoc.categorical
        decoded = Array{Union{decoded_eltype,Missing}}(missing, decoded_size...)
    else
        decoded = Array{Union{Tuple{decoded_eltype,decoded_eltype},Missing}}(missing, decoded_size...)
    end
    for ind in CartesianIndices(decoded)
        vec = encoded[ind,:]
        if !all(ismissing.(vec))
            decoded[ind] = decode(ohoc, vec)
        end
    end
    return decoded
end
