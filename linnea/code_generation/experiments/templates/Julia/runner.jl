using Test
using Logging
using MatrixGenerator

using LinearAlgebra.BLAS
BLAS.set_num_threads(1)

include("operand_generator.jl")
{0}
include("reference/naive.jl")
include("reference/recommended.jl")

matrices = operand_generator()

@info("Performing Test run...")
result_naive = collect(naive(map(copy, matrices)...))
result_recommended = collect(recommended(map(copy, matrices)...))
{1}
@test isapprox(result_recommended, result_naive)
{2}
@info("Test run performed successfully")


@info("Running Benchmarks...")
plotter = Benchmarker.Plot("julia_results_{3}.txt", ["algorithm"]);
{4}
Benchmarker.add_data(plotter, ["naive_julia"], Benchmarker.measure(20, naive, matrices...) );
Benchmarker.add_data(plotter, ["recommended_julia"], Benchmarker.measure(20, recommended, matrices...) );
Benchmarker.finish(plotter);
@info("Benchmarks complete")