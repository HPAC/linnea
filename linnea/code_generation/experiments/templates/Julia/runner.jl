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
result_naive = naive(matrices...)
result_recommended = recommended(matrices...)
@test isapprox(collect(result_recommended), collect(result_naive))
{1}
@info("Test run performed successfully")


@info("Running Benchmarks...")
plotter = Benchmarker.Plot("julia_results_{3}.txt", ["algorithm"]);
{2}
Benchmarker.add_data(plotter, ["naive_julia"], Benchmarker.measure(20, naive, matrices...) );
Benchmarker.add_data(plotter, ["recommended_julia"], Benchmarker.measure(20, recommended, matrices...) );
Benchmarker.finish(plotter);
@info("Benchmarks complete")