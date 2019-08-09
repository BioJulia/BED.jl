module BED

using Requires

function __init__()
    @require GenomicFeatures="899a7d2d-5c61-547b-bef9-6698a8d05446" include(joinpath(@__DIR__, "integrations", "GenomicFeatures.jl"))
end

import Automa
import Automa.RegExp: @re_str
import BGZFStreams
import BioCore
import BufferedStreams
import ColorTypes
import FixedPointNumbers: N0f8

include("record.jl")
include("reader.jl")
include("writer.jl")

end # module
