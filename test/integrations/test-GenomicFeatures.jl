using GenomicFeatures

# Testing strategy: there are two entirely separate intersection
# algorithms for IntervalCollection and IntervalStream. Here we test
# them both by checking that they agree by generating and intersecting
# random BED files.

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
            println(out, interval.seqname, "\t", interval.first - 1,
                    "\t", interval.last, "\t", interval.metadata, "\t",
                    1000, "\t", interval.strand)
        end
    end

end

n = 10000
Random.seed!(1234)
intervals_a = random_intervals(["one", "two", "three", "four", "five"],
                               1000000, n)
intervals_b = random_intervals(["one", "two", "three", "four", "five"],
                               1000000, n)

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
