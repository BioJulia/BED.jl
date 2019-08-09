using Test
using BED


using Distributions
import Random
import YAML
import ColorTypes: RGB
# import BioSequences: @dna_str, FASTA
import FixedPointNumbers: N0f8
import BGZFStreams

function get_bio_fmt_specimens(commit="222f58c8ef3e3480f26515d99d3784b8cfcca046")
    path = joinpath(dirname(@__FILE__), "BioFmtSpecimens")
    if !isdir(path)
        run(`git clone https://github.com/BioJulia/BioFmtSpecimens.git $(path)`)
    end
    cd(path) do
        #run(`git checkout $(commit)`)
    end
    return path
end

# Generate an array of n random Interval{Int} object. With sequence names
# samples from seqnames, and intervals drawn to lie in [1, maxpos].
function random_intervals(seqnames, maxpos::Int, n::Int)
    seq_dist = Categorical(length(seqnames))
    strand_dist = Categorical(2)
    length_dist = Normal(1000, 1000)
    intervals = Vector{Interval{Int}}(undef, n)
    for i in 1:n
        intlen = maxpos
        while intlen >= maxpos || intlen <= 0
            intlen = ceil(Int, rand(length_dist))
        end
        first = rand(1:maxpos-intlen)
        last = first + intlen - 1
        strand = rand(strand_dist) == 1 ? STRAND_POS : STRAND_NEG
        intervals[i] = Interval{Int}(seqnames[rand(seq_dist)],
                                     first, last, strand, i)
    end
    return intervals
end

@testset "BED" begin
    @testset "Record" begin
        record = BED.Record("chr1\t17368\t17436")
        @test BED.chrom(record) == "chr1"
        @test BED.chromstart(record) === 17369
        @test BED.chromend(record) === 17436
        @test !BED.hasname(record)
        @test !BED.hasscore(record)
        @test !BED.hasstrand(record)
        @test !BED.hasthickstart(record)
        @test !BED.hasthickend(record)
        @test !BED.hasitemrgb(record)
        @test !BED.hasblockcount(record)
        @test !BED.hasblocksizes(record)
        @test !BED.hasblockstarts(record)
        @test startswith(repr(record), "BED.Record:\n")
        @test string(record) == "chr1\t17368\t17436"

        record = BED.Record("chrXIII\t854794\t855293\tYMR292W\t0\t+\t854794\t855293\t0\t2\t22,395,\t0,104,")
        @test BED.chrom(record) == "chrXIII"
        @test BED.chromstart(record) === 854795
        @test BED.chromend(record) === 855293
        @test BED.name(record) == "YMR292W"
        @test BED.score(record) === 0
        @test BED.strand(record) === STRAND_POS
        @test BED.thickstart(record) === 854795
        @test BED.thickend(record) === 855293
        @test BED.itemrgb(record) === RGB(0, 0, 0)
        @test BED.blockcount(record) === 2
        @test BED.blocksizes(record) == [22, 395]
        @test BED.blockstarts(record) == [1, 105]
        @test startswith(repr(record), "BED.Record:\n")
        @test string(record) == "chrXIII\t854794\t855293\tYMR292W\t0\t+\t854794\t855293\t0\t2\t22,395,\t0,104,"

        record = BED.Record("chrX\t151080532\t151081699\tCHOCOLATE1\t0\t-\t151080532\t151081699\t255,127,36")
        @test BED.chrom(record) == "chrX"
        @test BED.chromstart(record) === 151080533
        @test BED.chromend(record) === 151081699
        @test BED.name(record) == "CHOCOLATE1"
        @test BED.score(record) === 0
        @test BED.strand(record) === STRAND_NEG
        @test BED.thickstart(record) === 151080533
        @test BED.thickend(record) === 151081699
        @test BED.itemrgb(record) === RGB(map(x -> reinterpret(N0f8, UInt8(x)), (255, 127, 36))...)
        @test !BED.hasblockcount(record)
        @test !BED.hasblocksizes(record)
        @test !BED.hasblockstarts(record)
        @test startswith(repr(record), "BED.Record:\n")
        @test string(record) == "chrX\t151080532\t151081699\tCHOCOLATE1\t0\t-\t151080532\t151081699\t255,127,36"
    end

    function check_bed_parse(filename)
        # Reading from a stream
        for interval in BED.Reader(open(filename))
        end

        # Reading from a regular file
        for interval in open(BED.Reader, filename)
        end

        # in-place parsing
        stream = open(BED.Reader, filename)
        entry = eltype(stream)()
        while !eof(stream)
            read!(stream, entry)
        end
        close(stream)

        # Check round trip
        output = IOBuffer()
        writer = BED.Writer(output)
        expected_entries = BED.Record[]
        for interval in open(BED.Reader, filename)
            write(writer, interval)
            push!(expected_entries, interval)
        end
        flush(writer)

        seekstart(output)
        read_entries = BED.Record[]
        for interval in BED.Reader(output)
            push!(read_entries, interval)
        end

        return expected_entries == read_entries
    end

    path = joinpath(get_bio_fmt_specimens(), "BED")
    for specimen in YAML.load_file(joinpath(path, "index.yml"))
        if "gzip" ∈ split(get(specimen, "tags", ""))
            # skip compressed files
            continue
        end
        valid = get(specimen, "valid", true)
        if valid
            @test check_bed_parse(joinpath(path, specimen["filename"]))
        else
            @test_throws Exception check_bed_parse(joinpath(path, specimen["filename"]))
        end
    end

    @testset "Integrations" begin
		include(joinpath("Integrations", "test-GenomicFeatures.jl"))
	end # testset Integrations

end
