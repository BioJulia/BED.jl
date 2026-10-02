# Benchmarks

These benchmarks compare BED.jl against [htslib](https://github.com/samtools/htslib), the C library behind `tabix` and `bgzip`.
htslib is loaded from [htslib_jll](https://github.com/JuliaBinaryWrappers/htslib_jll.jl) and called through the minimal `ccall` bindings in [`htslib.jl`](htslib.jl).

## Running

From the repository root, instantiate the environment once and then run the comparison:

```bash
julia --project=benchmark -e 'using Pkg; Pkg.instantiate()'
```

```bash
julia --project=benchmark benchmark/run.jl
```

The environment finds BED.jl in the parent directory through `[sources]`, which needs Julia 1.11 or later.
On older versions, run `Pkg.develop(path=".")` in the benchmark environment first.

The suite in [`benchmarks.jl`](benchmarks.jl) defines `SUITE` and also works with [PkgBenchmark](https://github.com/JuliaCI/PkgBenchmark.jl).
Set the environment variable `BED_BENCHMARK_RECORDS` to change the size of the data set (default 500000).

## Workloads

Synthetic BED6 and BED12 files are generated in a temporary directory, compressed with htslib's BGZF writer, and indexed with `tbx_index_build`.

- **scan**: read every record of the plain and the BGZF-compressed file.
  htslib reads with `tbx_readrec`, which splits each line and parses the sequence name, begin and end.
  BED.jl parses every field of each record, either by iterating over a `BED.Reader` (a new `Record` per line) or by reusing one `Record` with `read!`.
- **query**: find the records overlapping 100 random regions of 100 kb each with the tabix index.
  htslib uses `hts_itr_querys` and `hts_itr_next`, and BED.jl uses `eachoverlap`.
  The file and the index are opened once, outside the timed region.

Every variant returns the number of records it read and the sum of their interval lengths.
Each BED.jl variant is checked against htslib before it is benchmarked, and one that disagrees or throws is reported with a warning and left out of the suite.

## Caveats

htslib_jll provides htslib 1.19.1, which may lag behind the latest htslib release.
htslib does less work per record than BED.jl in the scan benchmarks, because it parses only the first three columns.
