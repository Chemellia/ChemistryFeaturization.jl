using ..ChemistryFeaturization.AbstractType: AbstractCodec
import ..ChemistryFeaturization.encode
import ..ChemistryFeaturization.decode
using Flux: onecold

export encode, decode, output_shape

"""
    OneHotOneCold(encode_f, decode_f, nbins, logspaced)

AbstractCodec type which uses a dummy variable (as defined in statistical literature), i.e., which employs
one-hot encoding and a one-cold decoding scheme.
"""
struct OneHotOneCold <: AbstractCodec
    categorical::Bool
    bins::Vector
end

output_shape(ohoc::OneHotOneCold) = ohoc.categorical ? length(ohoc.bins) : length(ohoc.bins) - 1

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

encode(ohoc::OneHotOneCold, vals::Vector) = encode.(Ref(ohoc), vals)

function decode(ohoc::OneHotOneCold, encoded)
    @assert size(encoded)[1] == output_shape(ohoc)
    local decoded
    if ohoc.categorical # return value
        decoded = onecold(encoded, ohoc.bins)
    else # return bounds (TODO: should this be a tuple or a vector..? I like tuple for distinguishing from encoded vectors, but it doesn't play so nice with broadcasting...)
        decoded = (onecold(encoded, ohoc.bins[1:end-1]), onecold(encoded, ohoc.bins[2:end]))
    end
    return decoded
end
