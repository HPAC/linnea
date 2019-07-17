using Test
using Logging
using MatrixGenerator

using LinearAlgebra.BLAS
BLAS.set_num_threads({num_threads})

include("operand_generator.jl")
{includes}
include("reference/naive.jl")
include("reference/recommended.jl")

function main()
    matrices = operand_generator()

    @info("Performing Test run...")
    result_naive = collect(naive(map(copy, matrices)...)[1])
    result_recommended = collect(recommended(map(copy, matrices)...)[1])
{tmp_results}
    @test isapprox(result_recommended, result_naive, rtol=1e-3)
{tests}
    @info("Test run performed successfully")


    @info("Running Benchmarks...")
    plotter = Benchmarker.Plot("julia_results_{experiment_name}", ["algorithm"; "threads"]);
{measurements}
    Benchmarker.add_data(plotter, ["naive_julia"; {num_threads}], Benchmarker.measure(20, naive, matrices...) );
    Benchmarker.add_data(plotter, ["recommended_julia"; {num_threads}], Benchmarker.measure(20, recommended, matrices...) );
    Benchmarker.finish(plotter);
    @info("Benchmarks complete")
end

main()
