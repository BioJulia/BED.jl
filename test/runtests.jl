using Test
using BED


using Distributions
using Documenter
using FormatSpecimens
using GenomicFeatures

import Random
import ColorTypes: RGB
import FixedPointNumbers: N0f8
import BGZFStreams

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

    dir_bed = path_of_format("BED")

    for specimen in list_valid_specimens("BED")

        if hastag(specimen, "gzip")
            # skip compressed files
            continue
        end

        filepath = joinpath(dir_bed, filename(specimen))

        @test check_bed_parse(filepath)

    end

    for specimen in list_invalid_specimens("BED")

        if hastag(specimen, "gzip")
            # skip compressed files
            continue
        end

        filepath = joinpath(dir_bed, filename(specimen))

        @test_throws Exception check_bed_parse(filepath)
    end

    # https://github.com/BioJulia/BED.jl/issues/13
    str = """

    # comment
    test1\t1\t9


    test2\t1\t9
    test3\t1\t9
        \t \t
    # comment
    # comment
    test4\t1\t9

    """

    @test collect(BED.Reader(IOBuffer(str))) == BED.Record.(string.("test", 1:4, "\t1\t9"))

    str = """
    # HOMER Peaks
    # Peak finding parameters:
    # tag directory = GM_tagdir
    #
    # total peaks = 158781
    # peak size = 153
    # peaks found using tags on both strands
    # minimum distance between peaks = 306
    # fragment length = 152
    # genome size = 2000000000
    # Total tags = 322230773.0
    # Total tags in peaks = 93989466.0
    # Approximate IP efficiency = 29.17%
    # tags per bp = 0.114786
    # expected tags per peak = 17.562
    # maximum tags considered per bp = 16.0
    # effective number of tags used for normalization = 10000000.0
    # Peaks have been centered at maximum tag pile-up
    # FDR rate threshold = 0.001000000
    # FDR effective poisson threshold = 1.591138e-05
    # FDR tag threshold = 38.0
    # number of putative peaks = 682969
    #
    # size of region used for local filtering = 10000
    # Fold over local region required = 4.00
    # Poisson p-value over local region required = 1.00e-04
    # Putative peaks filtered by local signal = 523066
    #
    # Maximum fold under expected unique positions for tags = 2.00
    # Putative peaks filtered for being too clonal = 1122
    #
    # cmd = findPeaks GM_tagdir -style factor -o GM.peaks.txt
    #
    # Column Headers:
    #PeakID\tchr\tstart\tend\tstrand\tNormalized Tag Count\tfocus ratio\tfindPeaks Score\tFold Change vs Local\tp-value vs Local\tClonal Fold Change
    chr21\t8401346\t8401499\tchr21-3\t1\t+
    chr21\t8445578\t8445731\tchr21-1\t1\t+
    chr21\t8218308\t8218461\tchr21-2\t1\t+
    """

    @test collect(BED.Reader(IOBuffer(str))) == [
        BED.Record("chr21\t8401346\t8401499\tchr21-3\t1\t+"),
        BED.Record("chr21\t8445578\t8445731\tchr21-1\t1\t+"),
        BED.Record("chr21\t8218308\t8218461\tchr21-2\t1\t+"),
    ]

    #=
    Testing strategy: there are two entirely separate intersection algorithms for IntervalCollection and IntervalStream.
    Here we test them both by checking that they agree by generating and intersecting random BED files.
    =#

    function check_intersection(filename_a, filename_b)
        ic_a = IntervalCollection{BED.Record}()
        open(BED.Reader, filename_a) do reader
            for record in reader
                push!(ic_a, Interval(record))
            end
        end

        ic_b = IntervalCollection{BED.Record}()
        open(BED.Reader, filename_b) do reader
            for record in reader
                push!(ic_b, Interval(record))
            end
        end

        # This is refactored out to close streams
        fa = open(BED.Reader, filename_a)
        fb = open(BED.Reader, filename_b)
        xs = sort(collect(eachoverlap(fa, fb)))
        close(fa)
        close(fb)

        ys = sort(collect(eachoverlap(ic_a, ic_b)))

        return xs == ys
    end

    function write_intervals(filename, intervals)
        open(filename, "w") do out
            for interval in sort(intervals)
                println(out, interval.seqname, "\t", interval.first - 1, "\t", interval.last, "\t", interval.metadata, "\t", 1000, "\t", interval.strand)
            end
        end

    end

    n = 10000
    Random.seed!(1234)
    intervals_a = random_intervals(["one", "two", "three", "four", "five"], 1000000, n)
    intervals_b = random_intervals(["one", "two", "three", "four", "five"], 1000000, n)

    mktempdir() do dir
        filename_a = joinpath(dir, "test_a.bed")
        filename_b = joinpath(dir, "test_b.bed")
        write_intervals(filename_a, intervals_a)
        write_intervals(filename_b, intervals_b)
        @test check_intersection(filename_a, filename_b)
    end

    @testset "eachoverlap" begin
        path = joinpath(path_of_format("BED"), "ws245Genes.WBGene.bed.bgz")
        stream = BGZFStreams.BGZFStream(path)
        reader = BED.Reader(stream, index=string(path, ".tbi"))
        for (interval, n_records) in [
                (Interval("chrII", 500_000:10_000_000), 38),
                (Interval("chrII", 1:70_000), 0),
                (Interval("chrIV", 500_000:10_000_000), 43),
                (Interval("chrI",  1_000_000:12_000_000), 27),]
            n = 0
            for record in eachoverlap(reader, interval)
                n += 1
            end
            @test n == n_records
        end
        @test isa(BED.Reader(path), BED.Reader)
    end

    # Include doctests.
    DocMeta.setdocmeta!(BED, :DocTestSetup, :(using BED); recursive=true)
    doctest(BED; manual = false)

end
