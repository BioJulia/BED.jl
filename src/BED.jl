module BED

using Automa: Automa, @re_str, @mark, @markpos, @relpos, @abspos
using Automa: onenter!, onexit!
import BGZFStreams
import BioGenerics
import ColorTypes
import FixedPointNumbers: N0f8

using GenomicFeatures
using Indexes
using TranscodingStreams

include("record.jl")
include("reader.jl")
include("writer.jl")

end # module
