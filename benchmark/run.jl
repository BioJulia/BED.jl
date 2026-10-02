# Run the BED.jl versus htslib benchmarks and print a comparison table.
#
# Usage, from the repository root:
#
#     julia --project=benchmark benchmark/run.jl
#
# Each row reports the median time and allocations of every variant, and the ratio of each BED.jl variant's median time to htslib's.

using BenchmarkTools
using Printf

include(joinpath(@__DIR__, "benchmarks.jl"))

println("Records: ", RECORD_COUNT, "; queries: ", QUERY_COUNT, " regions of ", QUERY_WIDTH, " bases")
println(Sys.cpu_info()[1].model, "; Julia ", VERSION, "; htslib_jll ", pkgversion(HTSlib.htslib_jll))
println()

tune!(SUITE)
results = run(SUITE; verbose = false)

function leaf_groups(group::BenchmarkGroup, path = String[])
    if all(value -> value isa BenchmarkTools.Trial, values(group))
        return [(join(path, " / "), group)]
    end
    return reduce(vcat, (leaf_groups(group[key], [path; key]) for key in sort!(collect(keys(group)))); init = [])
end

@printf("%-24s %-16s %12s %12s %10s\n", "benchmark", "variant", "median", "memory", "vs htslib")
for (name, group) in leaf_groups(results)
    reference = median(group["htslib"])
    for variant in sort!(collect(keys(group)); by = variant -> (variant == "htslib", variant))
        estimate = median(group[variant])
        ratio = variant == "htslib" ? "" : @sprintf("%.2f×", time(estimate) / time(reference))
        @printf("%-24s %-16s %12s %12s %10s\n", name, variant, BenchmarkTools.prettytime(time(estimate)), BenchmarkTools.prettymemory(memory(estimate)), ratio)
    end
end
