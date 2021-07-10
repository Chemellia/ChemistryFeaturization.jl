module Codec

include("simplecodec.jl")
export SimpleCodec

include("onehotonecold.jl")
export OneHotOneCold, build_onehot_vec

end
