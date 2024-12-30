using BenchmarkTools
push!(LOAD_PATH, "src")
using MyModule  # Import your module

suite = BenchmarkGroup()
suite["ec"] = BenchmarkGroup(["tag1","tag2"])
suite["ec"][10] = @benchmarkable expensive_computation(10)
suite["ec"][20] = @benchmarkable expensive_computation(20)
tune!(suite)
results = run(suite,verbose = true)
# Save the benchmark results to a JSON file
BenchmarkTools.save("benchmarks/benchmark_results.json", median(results))