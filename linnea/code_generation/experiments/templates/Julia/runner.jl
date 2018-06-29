using Base.Test
using MatrixGenerator

include("operand_generator.jl")
{0}
include("reference/naive.jl")
include("reference/recommended.jl")

matrices = operand_generator()
#test run
result_naive = naive(matrices...)
result_recommended = recommended(matrices...)
#@test result_recommended ≈ result_naive
#print("Norm: "*string(norm(result_naive - result_recommended))*"\\n")
#print("Naive(0, 0): "*string(result_naive[1, 1])*"\\n")
{1}

plotter = Benchmarker.Plotter.Plot{{Float64}}("julia_results_{3}.txt", ["algorithm"]);
{2}
Benchmarker.Plotter.add_data(plotter, ["naive_julia"], Benchmarker.measure(10, naive, matrices...) );
Benchmarker.Plotter.add_data(plotter, ["recommended_julia"], Benchmarker.measure(10, recommended, matrices...) );
Benchmarker.Plotter.finish(plotter);