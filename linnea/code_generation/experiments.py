
from ..algebra.properties import Property as properties
from .. import config
import textwrap
import os.path


op_gen_file_template = {config.Language.Julia: textwrap.dedent(
                            """
                            using MatrixGenerator

                            function operand_generator()
                            {}
                                return ({})
                            end
                            """
                        ),
                        config.Language.Matlab: textwrap.dedent(
                            """
                            function [out] = operand_generator()
                                import MatrixGenerator.*;
                            {}
                            end
                            """
                        ),
                        config.Language.Cpp: textwrap.dedent(
                            """
                            #include <generator/generator.hpp>

                            template<typename Gen>
                            decltype(auto) operand_generator(Gen && gen)
                            {{
                            {}
                                return std::make_tuple({});
                            }}
                            """
                        )
                        }
op_gen_line_template = {config.Language.Julia : "{name} = generate(({size}), [{properties}])",
                        config.Language.Cpp: "auto {name} = gen.generate({{{size}}}, {properties});",
                        config.Language.Matlab: "{name} = generate([{size}], {properties});"}

experiment_template = textwrap.dedent(
                            """
                            """
                            )

filename_extension = {config.Language.Julia: ".jl",
                      config.Language.Cpp: ".hpp",
                      config.Language.Matlab: ".m"
                      }

filename_extension_exec = {config.Language.Julia: ".jl",
                           config.Language.Cpp: ".cpp",
                           config.Language.Matlab: ".m"
                           }

operands_mapping_julia = {properties.SYMMETRIC: "Shape.Symmetric",
                          properties.DIAGONAL: "Shape.Diagonal",
                          properties.LOWER_TRIANGULAR: "Shape.LowerTriangular",
                          properties.UPPER_TRIANGULAR: "Shape.UpperTriangular",
                          properties.SPD: "Properties.SPD"
                         }

operands_mapping_matlab = {properties.SYMMETRIC: "Shape.Symmetric()",
                          properties.DIAGONAL: "Shape.Diagonal()",
                          properties.LOWER_TRIANGULAR: "Shape.LowerTriangular()",
                          properties.UPPER_TRIANGULAR: "Shape.UpperTriangular()",
                          properties.SPD: "Properties.SPD()"
                         }

operands_mapping_cpp = {properties.SYMMETRIC: "generator::shape::self_adjoint{}",
                        properties.DIAGONAL: "generator::shape::diagonal{}",
                        properties.LOWER_TRIANGULAR: "generator::shape::lower_triangular{}",
                        properties.UPPER_TRIANGULAR: "generator::shape::upper_triangular{}",
                        properties.SPD: "generator::property::spd{}"
                        }

operands_mapping = {config.Language.Julia: operands_mapping_julia,
                    config.Language.Cpp : operands_mapping_cpp,
                    config.Language.Matlab: operands_mapping_matlab}

random_operand =  {config.Language.Julia: ["Properties.Random", "Shape.General"],
                   config.Language.Cpp : ["generator::property::random{}"],
                   config.Language.Matlab: ["Properties.Random()", "Shape.General()"]}


benchmarker_code = {config.Language.Julia: textwrap.dedent(
                    """
                    using Base.Test
                    using MatrixGenerator
                    
                    include("operand_generator.jl")
                    {0}
                    include("algorithms/naive.jl")
                    include("algorithms/recommended.jl")
                    
                    matrices = operand_generator()
                    #test run
                    result_naive = naive(matrices...)
                    result_recommended = recommended(matrices...)
                    @test result_recommended ≈ result_naive
                    print("Norm: "*string(norm(result_naive - result_recommended))*"\\n")
                    print("Naive(0, 0): "*string(result_naive[1, 1])*"\\n")
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
                    function runner()
                        import MatrixGenerator.*;
                        % hardcoded directory with naive and recommended implementation
                        algorithms_dir = fullfile(fileparts(mfilename('fullpath')), 'algorithms');
                        addpath(algorithms_dir);
                        matrices = operand_generator();
                        naive_ = @() naive(matrices{{:}});
                        recommended_ = @() recommended(matrices{{:}});
                        naive_mat = naive_();
                        recomd_mat = recommended_();
                        fprintf('Norm(naive - recomd) %f\\n', norm(naive_mat - recomd_mat));
                        fprintf('Naive(0, 0): %f\\n', naive_mat(1, 1));
                        fprintf('Naive(0, 0): %f\\n', recomd_mat(1, 1));
                    
                        benchmarker = Benchmarker();
                        benchmarker.benchmark('naive_matlab', 20, naive_);
                        benchmarker.benchmark('recommended_matlab', 20, recommended_);
                        benchmarker.save('results.txt');
                    end
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
                        
                        #include "algorithms/naive_armadillo.hpp"
                        #include "algorithms/recommended_armadillo.hpp"
                        #include "algorithms/naive_eigen.hpp"
                        #include "algorithms/recommended_eigen.hpp"
                        #include "algorithms/naive_blaze.hpp"
                        
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
                            std::cout << "Armadillo norm(naive-recom): " << arma::norm(res_naive - res_recomm) << std::endl;
                            std::cout << "Armadillo naive(0,0): " <<res_naive(0, 0) << std::endl;
                        
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
                            std::cout << "Eigen norm(naive-recom): " << (res_naive - res_recomm).norm() << std::endl;
                            std::cout << "Eigen naive(0,0): " << res_naive(0, 0) << std::endl;
                        
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
                            std::cout << "Blaze naive(0, 0): " << res_naive(0, 0) << std::endl;
                        
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

algorithm_inclusion = {config.Language.Julia: "include(\"algorithms/algorithm{0}.jl\")",
                    config.Language.Cpp: ""}

algorithm_test = {config.Language.Julia: "@test algorithm{0}(map(MatrixGenerator.unwrap, matrices)...) ≈ result_naive",
                    config.Language.Cpp: ""}

algorithm_plot = {config.Language.Julia: "Benchmarker.Plotter.add_data(plotter, [\"generated{0}\"], Benchmarker.measure(20, algorithm{0}, map(MatrixGenerator.unwrap, matrices)...) );",
                    config.Language.Cpp: ""}

def map_operand(language, property):
    return operands_mapping.get(language).get(property)


# FIXME: somehow refactor this. it's ugly
def operand_generator_to_file(output_name, operands, output_str, language = config.Language.Julia, name_addition = ""):
    additional_name = "_{}".format(name_addition) if name_addition else ""
    file_name = os.path.join(config.output_path, config.language.name, output_name,
                             "operand_generator{}{}".format(additional_name, filename_extension.get(language)))
    directory_name = os.path.dirname(file_name)
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)
    output_file = open(file_name, "wt")
    op_gen_lines = []
    counter = 1
    for operand in operands:
        if language != config.Language.Matlab:
            replacement = {"name": operand.name}
        else:
            replacement = {"name": "out{{ {0} }}".format(counter)}
        property_replacements = []
        for prop in [properties.SYMMETRIC, properties.DIAGONAL, properties.LOWER_TRIANGULAR,
                     properties.UPPER_TRIANGULAR, properties.SPD]:
            if prop in operand.properties:
                property_replacements.append(map_operand(language, prop))
        if not properties.SPD in operand.properties:
            property_replacements.extend(random_operand.get(language))
        replacement["size"] = ",".join(map(str, operand.size))

        replacement["properties"] = ", ".join(property_replacements)
        op_gen_lines.append(op_gen_line_template.get(language).format(**replacement))

        # print(operand)
        counter += 1
    op_gen_lines = "\n".join(op_gen_lines)
    op_gen_file = op_gen_file_template.get(language).format(textwrap.indent(op_gen_lines, "    "), output_str)
    output_file.write(op_gen_file)
    output_file.close()


def benchmarker_to_file(output_name, language, algorithms_count=0):
    file_name = os.path.join(config.output_path, config.language.name, output_name,
                             "runner{}".format(filename_extension_exec.get(language)))
    output_file = open(file_name, "wt")
    inclusions = []
    tests = []
    plots = []
    incl_format = algorithm_inclusion.get(language)
    test_format = algorithm_test.get(language)
    plot_format = algorithm_plot.get(language)
    for i in range(algorithms_count):
        inclusions.append(incl_format.format(i))
        tests.append(test_format.format(i))
        plots.append(plot_format.format(i))
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

cmake_script = """
                cmake_minimum_required(VERSION 3.0)
                project(cpp_runner)
                
                find_package(LibTests REQUIRED)
                
                set(sources {0})
                
                foreach(source ${{sources}})
                    add_executable(${{source}} ${{source}}/runner.cpp)
                    target_link_libraries(${{source}} PRIVATE libtests)
                    set_target_properties(${{source}} PROPERTIES CXX_STANDARD 14)
                endforeach()
                """

def generate_cmake_script(paths):
    file_name = os.path.join(config.output_path, config.language.name, "CMakeLists.txt")
    output_file = open(file_name, "wt")
    names = " ".join(paths)
    output_file.write(cmake_script.format(names))
    output_file.close()
