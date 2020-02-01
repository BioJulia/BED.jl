# BED.jl
BED is a text-based file format for representing genomic annotations like genes, transcripts, and so on.
A BED file has tab-delimited and variable-length fields; the first three fields denoting a genomic interval are mandatory.

This is an example of RNA transcripts:
```
chr9	68331023	68424451	NM_015110	0	+
chr9	68456943	68486659	NM_001206	0	-
```

I/O tools for BED are provided by the `BED` module, which exports following three types:
* Reader type: `BED.Reader`
* Writer type: `BED.Writer`
* Element type: `BED.Record`

## Installation
You can install BED from the Julia REPL:

```julia
julia> Pkg.add("BED")
```

If you are interested in the cutting edge of the development, please check out the master branch to try new features before release.

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
