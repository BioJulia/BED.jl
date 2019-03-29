var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#BED.jl-1",
    "page": "Home",
    "title": "BED.jl",
    "category": "section",
    "text": "BED is a text-based file format for representing genomic annotations like genes, transcripts, and so on. A BED file has tab-delimited and variable-length fields; the first three fields denoting a genomic interval are mandatory.This is an example of RNA transcripts:chr9	68331023	68424451	NM_015110	0	+\nchr9	68456943	68486659	NM_001206	0	-I/O tools for BED are provided by the BED module, which exports following three types:Reader type: BED.Reader\nWriter type: BED.Writer\nElement type: BED.Record"
},

{
    "location": "#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "You can install BED from the Julia REPL:julia> Pkg.add(\"BED\")If you are interested in the cutting edge of the development, please check out the master branch to try new features before release."
},

{
    "location": "#Examples-1",
    "page": "Home",
    "title": "Examples",
    "category": "section",
    "text": "Here is a common workflow to iterate over all records in a BED file:# Import the BED module.\nusing BED\n\n# Open a BED file.\nreader = open(BED.Reader, \"data.bed\")\n\n# Iterate over records.\nfor record in reader\n    # Do something on record (see Accessors section).\n    chrom = BED.chrom(record)\n    # ...\nend\n\n# Finally, close the reader.\nclose(reader)If you repeatedly access records within specific ranges, it would be more efficient to construct an IntervalCollection object from a BED reader:using BED\nusing GenomicFeatures\n\n# Create an interval collection in memory.\nicol = open(BED.Reader, \"data.bed\") do reader\n    IntervalCollection(reader)\nend\n\n# Query overlapping records.\nfor interval in eachoverlap(icol, Interval(\"chrX\", 40001, 51500))\n    # A record is stored in the metadata field of an interval.\n    record = metadata(interval)\n    # ...\nend"
},

{
    "location": "api/#",
    "page": "API",
    "title": "API",
    "category": "page",
    "text": ""
},

{
    "location": "api/#BED.Reader",
    "page": "API",
    "title": "BED.Reader",
    "category": "type",
    "text": "BED.Reader(input::IO; index=nothing)\nBED.Reader(input::AbstractString; index=:auto)\n\nCreate a data reader of the BED file format.\n\nThe first argument specifies the data source. When it is a filepath that ends with .bgz, it is considered to be block compression file format (BGZF) and the function will try to find a tabix index file (<filename>.tbi) and read it if any. See http://www.htslib.org/doc/tabix.html for bgzip and tabix tools.\n\nArguments\n\ninput: data source\nindex: path to a tabix file\n\n\n\n\n\n"
},

{
    "location": "api/#BED.Writer",
    "page": "API",
    "title": "BED.Writer",
    "category": "type",
    "text": "BED.Writer(output::IO)\n\nCreate a data writer of the BED file format.\n\nArguments:\n\noutput: data sink\n\n\n\n\n\n"
},

{
    "location": "api/#BED.Record",
    "page": "API",
    "title": "BED.Record",
    "category": "type",
    "text": "BED.Record()\n\nCreate an unfilled BED record.\n\n\n\n\n\nBED.Record(data::Vector{UInt8})\n\nCreate a BED record object from data.\n\nThis function verifies and indexes fields for accessors. Note that the ownership of data is transferred to a new record object.\n\n\n\n\n\nBED.Record(str::AbstractString)\n\nCreate a BED record object from str.\n\nThis function verifies and indexes fields for accessors.\n\n\n\n\n\n"
},

{
    "location": "api/#BED.chrom",
    "page": "API",
    "title": "BED.chrom",
    "category": "function",
    "text": "chrom(record::Record)::String\n\nGet the chromosome name of record.\n\n\n\n\n\n"
},

{
    "location": "api/#BED.chromstart",
    "page": "API",
    "title": "BED.chromstart",
    "category": "function",
    "text": "chromstart(record::Record)::Int\n\nGet the starting position of record.\n\nNote that the first base is numbered 1.\n\n\n\n\n\n"
},

{
    "location": "api/#BED.chromend",
    "page": "API",
    "title": "BED.chromend",
    "category": "function",
    "text": "chromend(record::Record)::Int\n\nGet the end position of record.\n\n\n\n\n\n"
},

{
    "location": "api/#BED.name",
    "page": "API",
    "title": "BED.name",
    "category": "function",
    "text": "name(record::Record)::String\n\nGet the name of record.\n\n\n\n\n\n"
},

{
    "location": "api/#BED.score",
    "page": "API",
    "title": "BED.score",
    "category": "function",
    "text": "score(record::Record)::Int\n\nGet the score between 0 and 1000.\n\n\n\n\n\n"
},

{
    "location": "api/#BED.strand",
    "page": "API",
    "title": "BED.strand",
    "category": "function",
    "text": "strand(record::Record)::GenomicFeatures.Strand\n\nGet the strand of record.\n\n\n\n\n\n"
},

{
    "location": "api/#BED.thickstart",
    "page": "API",
    "title": "BED.thickstart",
    "category": "function",
    "text": "thickstart(record::Record)::Int\n\nGet the starting position at which record is drawn thickly.\n\nNote that the first base is numbered 1.\n\n\n\n\n\n"
},

{
    "location": "api/#BED.thickend",
    "page": "API",
    "title": "BED.thickend",
    "category": "function",
    "text": "thickend(record::Record)::Int\n\nGet the end position at which record is drawn thickly.\n\n\n\n\n\n"
},

{
    "location": "api/#BED.itemrgb",
    "page": "API",
    "title": "BED.itemrgb",
    "category": "function",
    "text": "itemrgb(record::Record)::ColorTypes.RGB\n\nGet the RGB value of record.\n\nThe return type is defined in ColorTypes.jl.\n\n\n\n\n\n"
},

{
    "location": "api/#BED.blockcount",
    "page": "API",
    "title": "BED.blockcount",
    "category": "function",
    "text": "blockcount(record::Record)::Int\n\nGet the number of blocks (exons) in record.\n\n\n\n\n\n"
},

{
    "location": "api/#BED.blocksizes",
    "page": "API",
    "title": "BED.blocksizes",
    "category": "function",
    "text": "blocksizes(record::Record)::Vector{Int}\n\nGet the block (exon) sizes of record.\n\n\n\n\n\n"
},

{
    "location": "api/#BED.blockstarts",
    "page": "API",
    "title": "BED.blockstarts",
    "category": "function",
    "text": "blockstarts(record::Record)::Vector{Int}\n\nGet the block (exon) starts of record.\n\nNote that the first base is numbered 1.\n\n\n\n\n\n"
},

{
    "location": "api/#API-1",
    "page": "API",
    "title": "API",
    "category": "section",
    "text": "BED.Reader\nBED.Writer\nBED.Record\nBED.chrom\nBED.chromstart\nBED.chromend\nBED.name\nBED.score\nBED.strand\nBED.thickstart\nBED.thickend\nBED.itemrgb\nBED.blockcount\nBED.blocksizes\nBED.blockstarts"
},

]}
