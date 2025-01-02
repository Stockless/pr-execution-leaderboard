using BenchmarkTools
using MyPkg

# Get GitHub username from environment variable
github_username = get(ENV, "GITHUB_ACTOR", "UNKNOWN")

suite = BenchmarkGroup()
suite["function1"][github_username] = @benchmarkable expensive_computation(30000)
tune!(suite)
results = run(suite, verbose=true)

# Save the benchmark results to a JSON file
BenchmarkTools.save("benchmarks/benchmark_results.json", median(results), std(results))