using ..ChemistryFeaturization.AbstractType: AbstractCodec

"""
    OneHotOneCold(encode_f, decode_f, nbins, logspaced)

AbstractCodec type which uses a dummy variable (as defined in statistical literature), i.e., which employs
one-hot encoding and a one-cold decoding scheme.
"""
struct OneHotOneCold <: AbstractCodec
    encode_f::Function
    decode_f::Function
    nbins::Integer
    logspaced::Bool
end

"A flexible version of Flux.onehot that can handle both categorical and continuous-valued encoding."
function build_onehot_vec(val, bins, categorical)
    local bin_index , onehot_vec
    if categorical
        onehot_vec = [0.0 for i = 1:length(bins)]
        bin_index = findfirst(isequal(val), bins)
    else
        @assert eltype(bins) <: Number "Your bins aren't numbers...are you sure you didn't mean for this feature to be categorical?"
        @assert bins[1] <= val <= bins[end] "The value $val is outside the range of bins $bins"
        onehot_vec = [0.0 for i = 1:(length(bins)-1)]
        bin_index = searchsorted(bins, val).stop
        if bin_index == length(bins) # got the max value
            bin_index = bin_index - 1
        elseif isapprox(val, bins[1]) # sometimes will get 0 if this doesn't get checked
            bin_index = 1
        end
    end
    onehot_vec[bin_index] = 1.0
    return onehot_vec
end