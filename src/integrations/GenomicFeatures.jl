function GenomicFeatures.Interval(record::Record)
    name = BioCore.seqname(record)
    lpos = BioCore.leftposition(record)
    rpos = BioCore.rightposition(record)
    strd = hasstrand(record) ? GenomicFeatures.strand(record) : GenomicFeatures.STRAND_BOTH
    return GenomicFeatures.Interval(name, lpos, rpos, strd, record)
end

function Base.convert(::Type{GenomicFeatures.Interval}, record::Record)
    return GenomicFeatures.Interval(record)
end

function Base.convert(::Type{GenomicFeatures.Interval{Record}}, record::Record)
    return convert(GenomicFeatures.Interval, record)
end

function GenomicFeatures.eachoverlap(reader::Reader, interval::GenomicFeatures.Interval)
    if reader.index === nothing
        throw(ArgumentError("index is null"))
    end
    return GenomicFeatures.Indexes.TabixOverlapIterator(reader, interval)
end

function GenomicFeatures.strand(record::Record)
    return strand(record)
end
