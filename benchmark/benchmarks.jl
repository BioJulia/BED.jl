# BED.jl benchmarks against htslib
# ================================
#
# Defines `SUITE`, a BenchmarkTools.BenchmarkGroup compatible with PkgBenchmark.
# Each benchmark has "BED.jl" and "htslib" variants doing equivalent work, and every variant returns a checksum (record count and summed interval length).
# When this file is loaded, each BED.jl variant is checked against htslib's checksum, and one that disagrees or throws is reported and left out.
#
# The number of records in the synthetic data set is set by the environment variable `BED_BENCHMARK_RECORDS` (default 500_000).

using BED
using BenchmarkTools
using GenomicFeatures
using Random

include(joinpath(@__DIR__, "htslib.jl"))
using .HTSlib

const RECORD_COUNT = parse(Int, get(ENV, "BED_BENCHMARK_RECORDS", "500000"))
const QUERY_COUNT = 100
const QUERY_WIDTH = 100_000
const CHROMOSOMES = ["chr1", "chr2", "chr3", "chr4", "chr5"]
const STRANDS = ('+', '-', '.')

"""
    write_bed(path, record_count; columns, rng)

Write `record_count` random BED records sorted by chromosome and start position, as tabix requires.
There are no header or comment lines, which `tbx_readrec` would reject.
"""
function write_bed(path::AbstractString, record_count::Integer; columns::Integer = 6, rng::AbstractRNG = Xoshiro(1234))
    records_per_chromosome = cld(record_count, length(CHROMOSOMES))
    open(path, "w") do io
        written = 0
        for chromosome in CHROMOSOMES
            position = 0
            for _ in 1:min(records_per_chromosome, record_count - written)
                position += rand(rng, 0:200)
                width = rand(rng, 1:2_000)
                print(io, chromosome, '\t', position, '\t', position + width)
                if columns >= 6
                    print(io, "\tfeature", written + 1, '\t', rand(rng, 0:1000), '\t', rand(rng, STRANDS))
                end
                if columns >= 12
                    block_size = cld(width, 2)
                    print(io, '\t', position, '\t', position + width, "\t255,0,0\t2\t", block_size, ',', width - block_size, ",\t0,", block_size, ',')
                end
                println(io)
                written += 1
            end
        end
    end
    return path
end

"""
    random_regions(record_count; rng)

Draw query regions of `QUERY_WIDTH` bases spread over the span of the synthetic chromosomes.
"""
function random_regions(record_count::Integer; rng::AbstractRNG = Xoshiro(5678))
    # Records start on average every 100 bases.
    span = 100 * cld(record_count, length(CHROMOSOMES))
    return map(1:QUERY_COUNT) do _
        chromosome = rand(rng, CHROMOSOMES)
        first = rand(rng, 1:max(1, span - QUERY_WIDTH))
        return GenomicFeatures.Interval(chromosome, first, first + QUERY_WIDTH - 1)
    end
end

region_string(interval::GenomicFeatures.Interval) = string(GenomicFeatures.seqname(interval), ':', first(interval), '-', last(interval))

# BED.jl workloads, mirroring HTSlib.scan and HTSlib.query.
# `BED.Reader(path)` forwards `index = :auto` to `BED.Reader(::IO)` for plain files, which rejects it, so plain files are opened here.
open_reader(path::AbstractString) = endswith(path, ".bgz") ? BED.Reader(path) : BED.Reader(open(path))

# BED.jl positions are 1-based and inclusive, so `chromend - chromstart + 1` matches htslib's `end - begin`.

function bed_scan_iterate(path::AbstractString)
    count = 0
    total_length = 0
    reader = open_reader(path)
    for record in reader
        count += 1
        total_length += BED.chromend(record) - BED.chromstart(record) + 1
    end
    close(reader)
    return count, total_length
end

function bed_scan_inplace(path::AbstractString)
    count = 0
    total_length = 0
    reader = open_reader(path)
    record = BED.Record()
    while BED.BioGenerics.IO.tryread!(reader, empty!(record)) !== nothing
        count += 1
        total_length += BED.chromend(record) - BED.chromstart(record) + 1
    end
    close(reader)
    return count, total_length
end

function bed_query(reader::BED.Reader, region::GenomicFeatures.Interval)
    count = 0
    total_length = 0
    for record in eachoverlap(reader, region)
        count += 1
        total_length += BED.chromend(record) - BED.chromstart(record) + 1
    end
    return count, total_length
end

function bed_queries(reader::BED.Reader, regions)
    count = 0
    total_length = 0
    for region in regions
        (region_count, region_length) = bed_query(reader, region)
        count += region_count
        total_length += region_length
    end
    return count, total_length
end

function htslib_queries(indexed::HTSlib.IndexedFile, regions)
    count = 0
    total_length = 0
    for region in regions
        (region_count, region_length) = HTSlib.query(indexed, region)
        count += region_count
        total_length += region_length
    end
    return count, total_length
end

# Data
# ----

const DATA_DIRECTORY = mktempdir()

const DATA = map((bed6 = 6, bed12 = 12)) do columns
    plain = write_bed(joinpath(DATA_DIRECTORY, "bed$(columns).bed"), RECORD_COUNT; columns)
    # BED.Reader recognises BGZF by the ".bgz" extension and finds the tabix index at "<path>.tbi".
    compressed = HTSlib.bgzip(plain, plain * ".bgz")
    HTSlib.tabix_index(compressed)
    # Open handles are shared by the benchmarks so that loading the index is not timed.
    reader = BED.Reader(compressed)
    indexed = HTSlib.IndexedFile(compressed)
    return (; plain, compressed, reader, indexed)
end
atexit(() -> foreach(files -> (close(files.reader); close(files.indexed)), DATA))

const REGIONS = random_regions(RECORD_COUNT)
const REGION_STRINGS = map(region_string, REGIONS)

"""
    agrees(description, workload, expected)

Run `workload()` and check that it returns the htslib checksum `expected`.
A BED.jl variant that errors or disagrees is reported and left out of the suite, so that a bug does not stop the remaining benchmarks.
"""
function agrees(description::AbstractString, workload, expected)
    actual = try
        workload()
    catch exception
        @warn "Skipping $(description): it threw an exception" exception
        return false
    end
    if actual != expected
        @warn "Skipping $(description): it returned (count, total length) $(actual) but htslib returned $(expected)"
        return false
    end
    return true
end

# Suite
# -----
#
# Both implementations are checked to do the same work before their benchmarks are added.

const SUITE = BenchmarkGroup()
SUITE["scan"] = BenchmarkGroup()
SUITE["query"] = BenchmarkGroup()

for (name, files) in pairs(DATA)
    SUITE["scan"][string(name)] = BenchmarkGroup()
    for (compression, path) in (("plain", files.plain), ("bgzf", files.compressed))
        group = SUITE["scan"][string(name)][compression] = BenchmarkGroup()
        expected = HTSlib.scan(path, files.indexed)
        expected[1] == RECORD_COUNT || error("htslib read $(expected[1]) of $(RECORD_COUNT) records from $(path)")
        group["htslib"] = @benchmarkable HTSlib.scan($path, $(files.indexed))
        if agrees("scan/$(name)/$(compression)/BED.jl iterate", () -> bed_scan_iterate(path), expected)
            group["BED.jl iterate"] = @benchmarkable bed_scan_iterate($path)
        end
        if agrees("scan/$(name)/$(compression)/BED.jl read!", () -> bed_scan_inplace(path), expected)
            group["BED.jl read!"] = @benchmarkable bed_scan_inplace($path)
        end
    end
    group = SUITE["query"][string(name)] = BenchmarkGroup()
    expected = htslib_queries(files.indexed, REGION_STRINGS)
    group["htslib"] = @benchmarkable htslib_queries($(files.indexed), $REGION_STRINGS)
    if agrees("query/$(name)/BED.jl", () -> bed_queries(files.reader, REGIONS), expected)
        group["BED.jl"] = @benchmarkable bed_queries($(files.reader), $REGIONS)
    end
end
