module Codec

@enum EncodeOrDecode ENCODE DECODE
export EncodeOrDecode

include("simplecodec.jl")
export SimpleCodec

include("onehotonecold.jl")
export OneHotOneCold

end
