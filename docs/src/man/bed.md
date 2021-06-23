# BED
BED is a text-based file format for representing genomic annotations like genes, transcripts, and so on.
A BED file has tab-delimited and variable-length fields; the first three fields denoting a genomic interval are mandatory.

This is an example of RNA transcripts:
```
chr9	68331023	68424451	NM_015110	0	+
chr9	68456943	68486659	NM_001206	0	-
```

The `BED` package supports I/O for BED by providing the following three types:
* Reader type: `BED.Reader`
* Writer type: `BED.Writer`
* Element type: `BED.Record`

## Examples

Here is a common workflow to iterate over all records in a BED file:
```julia
# Import the BED module.
using BED

# Open a BED file.
reader = open(BED.Reader, "data.bed")

# Iterate over records.
for record in reader
    # Do something on record (see Accessors section).
    chrom = BED.chrom(record)
    # ...
end

# Finally, close the reader.
close(reader)
```

The iterator interface demonstrated above allocates an object for each record and that may be a bottleneck of reading data from a file.
In-place reading reuses a pre-allocated object for every record and less memory allocation happens in reading:

```julia
# Import the BED module.
using BED

# Open a BED file.
reader = open(BED.Reader, "data.bed")

# Pre-allocate record.
record = BED.Record()
while !eof(reader)
    empty!(record)
    read!(reader, record)
    # do something
end

# Finally, close the reader.
close(reader)
```

If you repeatedly access records within specific ranges, it would be more efficient to construct an `IntervalCollection` object from a BED reader:
```julia
using BED
using GenomicFeatures

# Create an interval collection in memory.
icol = open(BED.Reader, "data.bed") do reader
    IntervalCollection(reader)
end

# Query overlapping records.
for interval in eachoverlap(icol, Interval("chrX", 40001, 51500))
    # A record is stored in the metadata field of an interval.
    record = metadata(interval)
    # ...
end
```
