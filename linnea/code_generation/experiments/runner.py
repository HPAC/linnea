
import textwrap
import os.path

from ... import config

filename_extension_exec = {config.Language.Julia: ".jl",
                           config.Language.Cpp: ".cpp",
                           config.Language.Matlab: ".m"
                           }

algorithm_inclusion = {config.Language.Julia: """include("{0}/{1}.jl")""",
                    config.Language.Cpp: ""}

algorithm_test = {config.Language.Julia: "#@test {0}(map(MatrixGenerator.unwrap, matrices)...) ≈ result_naive",
                    config.Language.Cpp: ""}

algorithm_plot = {config.Language.Julia: """Benchmarker.Plotter.add_data(plotter, ["{0}"], Benchmarker.measure(20, {0}, map(MatrixGenerator.unwrap, matrices)...) );""",
                    config.Language.Cpp: ""}

benchmarker_code = {config.Language.Julia: textwrap.dedent(
                    """
                    using Test
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
                    """
                    ),
                    config.Language.Matlab: textwrap.dedent(
                    """
                    import MatrixGenerator.*;
                    % hardcoded directory with naive and recommended implementation
                    algorithms_dir = fullfile(fileparts(mfilename('fullpath')), 'algorithms');
                    addpath(algorithms_dir);
                    matrices = operand_generator();
                    naive_ = @() naive(matrices{{:}});
                    recommended_ = @() recommended(matrices{{:}});
                    naive_mat = naive_();
                    recomd_mat = recommended_();
                    %fprintf('Norm(naive - recomd) %f\\n', norm(naive_mat - recomd_mat));
                    %fprintf('Naive(0, 0): %f\\n', naive_mat(1, 1));
                    %fprintf('Naive(0, 0): %f\\n', recomd_mat(1, 1));
                
                    benchmarker = Benchmarker();
                    benchmarker.benchmark('naive_matlab', 10, naive_);
                    benchmarker.benchmark('recommended_matlab', 10, recommended_);
                    benchmarker.save('matlab_results_{0}.txt');
                    """
                    ),
                    config.Language.Cpp: textwrap.dedent(
                        """
                        #include <array>
                        #include <string>
                        #include <utility>
                        #include <fstream>
                        
                        #include <benchmarker/benchmarker.hpp>
                        #include <libraries.hpp>
                        #include <generator/generator.hpp>
                        
                        #include "reference/naive_armadillo.hpp"
                        #include "reference/recommended_armadillo.hpp"
                        #include "reference/naive_eigen.hpp"
                        #include "reference/recommended_eigen.hpp"
                        #include "reference/naive_blaze.hpp"
                        
                        #include "operand_generator.hpp"
                        
                        typedef double float_type;
                        
                        template<typename Function, typename Tuple, std::size_t... I>
                        decltype(auto) call(Function && f, Tuple && t, std::index_sequence<I...>)
                        {{
                            return f(std::get<I>(t)...);
                        }}
                        
                        int main(int argc, char ** argv)
                        {{
                            std::cout << "Test runner for algorithm {0}\\n";
                            linalg_tests::basic_benchmarker<std::chrono::duration<float>> benchmark;
                            benchmark.set_cache_size(30 * 1024 * 1024);
                            {{
                            generator::generator<library::arma, float_type> arma_gen{{100}};
                            auto matrices = operand_generator(arma_gen);
                            constexpr std::size_t tuple_length = std::tuple_size<decltype(matrices)>::value;
                            typedef std::make_index_sequence<tuple_length> Indices;
                            auto res_naive = call(naive_armadillo{{}}, matrices, Indices{{}});
                            auto res_recomm = call(recommended_armadillo{{}}, matrices, Indices{{}});
                            //std::cout << "Armadillo norm(naive-recom): " << arma::norm(res_naive - res_recomm) << std::endl;
                            //std::cout << "Armadillo naive(0,0): " <<res_naive(0, 0) << std::endl;
                        
                            benchmark.run(20, [&]() {{
                                    call(naive_armadillo{{}}, matrices, Indices{{}});
                                    }});
                            benchmark.run(20, [&]() {{
                                    call(recommended_armadillo{{}}, matrices, Indices{{}});
                                    }});
                            }}
                            {{ 
                            generator::generator<library::eigen, float_type> eigen_gen{{100}};
                            auto matrices = operand_generator(eigen_gen);
                            constexpr std::size_t tuple_length = std::tuple_size<decltype(matrices)>::value;
                            typedef std::make_index_sequence<tuple_length> Indices;
                            auto res_naive = call(naive_eigen{{}}, matrices, Indices{{}});
                            auto res_recomm = call(recommended_eigen{{}}, matrices, Indices{{}});
                            //std::cout << "Eigen norm(naive-recom): " << (res_naive - res_recomm).norm() << std::endl;
                            //std::cout << "Eigen naive(0,0): " << res_naive(0, 0) << std::endl;
                        
                            benchmark.run(20, [&]() {{
                                    call(naive_eigen{{}}, matrices, Indices{{}});
                                    }});
                            benchmark.run(20, [&]() {{
                                    call(recommended_eigen{{}}, matrices, Indices{{}});
                                    }});
                            }}
                            {{
                            generator::generator<library::blaze, float_type> blaze_gen{{100}};
                            auto matrices = operand_generator(blaze_gen);
                            constexpr std::size_t tuple_length = std::tuple_size<decltype(matrices)>::value;
                            typedef std::make_index_sequence<tuple_length> Indices;
                            auto res_naive = call(naive_blaze{{}}, matrices, Indices{{}});
                            //std::cout << "Blaze naive(0, 0): " << res_naive(0, 0) << std::endl;
                        
                            benchmark.run(20, [&]() {{
                                    call(naive_blaze{{}}, matrices, Indices{{}});
                                    }});
                            }}
                        
                            std::array<std::string, 5> labels{{ {{"naive_armadillo", "recommended_armadillo",
                                "naive_eigen", "recommended_eigen", "naive_blaze"}} }};
                            std::ofstream file("cpp_results_{0}.txt");
                            benchmark.output_results(file, labels);
                            file.close();
                        }}
                        """
                    )
                }



def runner_to_file(output_name, language, algorithms=None):
    file_name = os.path.join(config.output_path, output_name, language.name,
                             "runner{}".format(filename_extension_exec.get(language)))
    output_file = open(file_name, "wt", encoding='utf-8')
    inclusions = []
    tests = []
    plots = []
    incl_format = algorithm_inclusion.get(language)
    test_format = algorithm_test.get(language)
    plot_format = algorithm_plot.get(language)
    if algorithms:
        for subdir_name, algorithm_name in algorithms:
            inclusions.append(incl_format.format(subdir_name, algorithm_name))
            tests.append(test_format.format(algorithm_name))
            plots.append(plot_format.format(algorithm_name))
    if language == config.Language.Julia:
        op_gen_file = benchmarker_code.get(language).format(
            "\n".join(inclusions),
            "\n".join(tests),
            "\n".join(plots),
            output_name
        )
    else:
        op_gen_file = benchmarker_code.get(language).format(output_name)
    output_file.write(op_gen_file)
    output_file.close()


cmake_script = textwrap.dedent(
                """
                cmake_minimum_required(VERSION 3.0)
                project({0})
                
                find_package(LibTests REQUIRED)
                
                add_executable({0} runner.cpp)
                target_link_libraries({0} PRIVATE libtests)
                set_target_properties({0} PROPERTIES CXX_STANDARD 14)
                """)


def generate_cmake_script(output_name):
    file_name = os.path.join(config.output_path, output_name, config.Language.Cpp.name, "CMakeLists.txt")
    output_file = open(file_name, "wt")
    output_file.write(cmake_script.format(output_name))
    output_file.close()


def generate_runner(output_name, algorithms):

    runner_to_file(output_name, algorithms=algorithms, language=config.Language.Julia)
    runner_to_file(output_name, language=config.Language.Matlab)
    runner_to_file(output_name, language=config.Language.Cpp)
    generate_cmake_script(output_name)
