module Codec

@enum EncodeOrDecode ENCODE DECODE
export EncodeOrDecode

include("OneHotOneCold.jl")
export OneHotOneCold, build_onehot_vec

end
