
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
                        config.Language.Cpp: "auto {name} = gen.generate({{{size}}}, {properties});"}

experiment_template = textwrap.dedent(
                            """
                            """
                            )

filename_extension = {config.Language.Julia : ".jl",
                      config.Language.Cpp: ".hpp"}

operands_mapping_julia = {properties.SYMMETRIC: "Shape.Symmetric",
                          properties.DIAGONAL: "Shape.Diagonal",
                          properties.LOWER_TRIANGULAR: "Shape.LowerTriangular",
                          properties.UPPER_TRIANGULAR: "Shape.UpperTriangular",
                          properties.SPD: "Properties.SPD"
                         }

operands_mapping_cpp = {properties.SYMMETRIC: "generator::shape::self_adjoint{}",
                        properties.DIAGONAL: "generator::shape::diagonal{}",
                        properties.LOWER_TRIANGULAR: "generator::shape::lower_triangular{}",
                        properties.UPPER_TRIANGULAR: "generator::shape::upper_triangular{}",
                        properties.SPD: "generator::property::spd{}"
                        }
operands_mapping = {config.Language.Julia: operands_mapping_julia,
                    config.Language.Cpp : operands_mapping_cpp}

random_operand =  {config.Language.Julia: ["Properties.Random", "Shape.General"],
                    config.Language.Cpp : ["generator::property::random{}"]}


benchmarker_code = {config.Language.Julia: textwrap.dedent(
                    """
                    using Base.Test
                    using MatrixGenerator
                    
                    include("operand_generator.jl")
                    have_algorithm = isfile("algorithms/algorithm0.jl")
                    if have_algorithm
                        include("algorithms/algorithm0.jl")
                    end
                    {0}
                    include("algorithms/naive.jl")
                    include("algorithms/recommended.jl")
                    
                    matrices = operand_generator()
                    #test run
                    result_naive = naive(matrices...)
                    result_recommended = recommended(matrices...)
                    @test result_recommended ≈ result_naive
                    {1}
                    
                    plotter = Benchmarker.Plotter.Plot{{Float64}}("julia_results.txt", ["algorithm"]);
                    {2}
                    Benchmarker.Plotter.add_data(plotter, ["naive"], Benchmarker.measure(10, naive, matrices...) );
                    Benchmarker.Plotter.add_data(plotter, ["recommended"], Benchmarker.measure(10, recommended, matrices...) );
                    Benchmarker.Plotter.finish(plotter);
                    """
                ),
                config.Language.Cpp: textwrap.dedent(
                    """
                    """
                )
                }

algorithm_inclusion = {config.Language.Julia: "include(\"algorithms/algorithm{0}.jl\")",
                    config.Language.Cpp: ""}

algorithm_test = {config.Language.Julia: "@test algorithm{0}(map(MatrixGenerator.unwrap, matrices)...) ≈ result_naive",
                    config.Language.Cpp: ""}

algorithm_plot = {config.Language.Julia: "Benchmarker.Plotter.add_data(plotter, [\"algorithm{0}\"], Benchmarker.measure(20, algorithm{0}, map(MatrixGenerator.unwrap, matrices)...) );",
                    config.Language.Cpp: ""}

def map_operand(language, property):
    return operands_mapping.get(language).get(property)

def operand_generator_to_file(output_name, operands, output_str, language = config.Language.Julia, name_addition = ""):
    additional_name = "_{}".format(name_addition) if name_addition else ""
    file_name = os.path.join(config.output_path, config.language.name, output_name,
                             "operand_generator{}{}".format(additional_name, filename_extension.get(language)))
    directory_name = os.path.dirname(file_name)
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)
    output_file = open(file_name, "wt")
    op_gen_lines = []
    for operand in operands:
        replacement = {"name": operand.name}
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
    op_gen_lines = "\n".join(op_gen_lines)
    # print(op_gen_lines)
    op_gen_file = op_gen_file_template.get(language).format(textwrap.indent(op_gen_lines, "    "), output_str)
    output_file.write(op_gen_file)
    output_file.close()

def benchmarker_to_file(output_name, algorithms_count, language):
    file_name = os.path.join(config.output_path, config.language.name, output_name,
                             "runner{}".format(filename_extension.get(language)))
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
    op_gen_file = benchmarker_code.get(language).format(
        "\n".join(inclusions),
        "\n".join(tests),
        "\n".join(plots)
    )
    output_file.write(op_gen_file)
    output_file.close()

