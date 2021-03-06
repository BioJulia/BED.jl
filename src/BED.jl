module BED

import Automa
import Automa.RegExp: @re_str
import BGZFStreams
import BioCore
import BufferedStreams
import ColorTypes
import FixedPointNumbers: N0f8

using GenomicFeatures
using Indexes

include("record.jl")
include("reader.jl")
include("writer.jl")

end # module
