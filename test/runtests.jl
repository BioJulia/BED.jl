using Test
using BED

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
        @test startswith(repr(record), "GenomicFeatures.BED.Record:\n")
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
        @test startswith(repr(record), "GenomicFeatures.BED.Record:\n")
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
        @test startswith(repr(record), "GenomicFeatures.BED.Record:\n")
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
        path = joinpath(get_bio_fmt_specimens(), "BED", "ws245Genes.WBGene.bed.bgz")
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
end
