module Codecs

@enum EncodeOrDecode ENCODE DECODE
export EncodeOrDecode

include("OneHotOneCold.jl")
export OneHotOneCold

end
