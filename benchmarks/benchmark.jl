using BenchmarkTools
push!(LOAD_PATH, "src")
using MyModule  # Import your module

# Benchmark the expensive_computation function
@benchmark expensive_computation(30)  # Benchmark for n = 30

# Benchmark the fast_function
@benchmark fast_function(100)
