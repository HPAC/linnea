
from ..algebra.properties import Property as properties
from .. import config
import textwrap
import os.path

algorithm_template = textwrap.dedent(
                        """
                        using Base.LinAlg.BLAS
                        using Base.LinAlg

                        function {}({})
                        {}
                            return ({})
                        end
                        """)

op_gen_file_template = textwrap.dedent(
                        """
                        push!(LOAD_PATH, "/Users/henrik/Promotion/code/linalg_tests/julia/")
                        using Benchmarker

                        function operand_generator()
                        {}
                            return ({})
                        end
                        """)

op_gen_line_template = """{name} = generate(Shape.{shape}({size}), Set([{properties}]))"""

experiment_template = textwrap.dedent(
                            """
                            """
                            )


def experiment_to_file(experiment_name, algorithm_name, algorithm, input, output):
    file_name = os.path.join(config.experiments_path, config.language.name, experiment_name, "algorithms", "{}.jl".format(algorithm_name))
    directory_name = os.path.dirname(file_name)
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)
    output_file = open(file_name, "wt")
    experiment_str = algorithm_template.format(algorithm_name, input, textwrap.indent(algorithm, "    "), output)
    output_file.write(experiment_str)
    output_file.close()
    # print(file_name)

def operand_generator_to_file(experiment_name, operands, output_str):
    file_name = os.path.join(config.experiments_path, config.language.name, experiment_name, "operand_generator.jl")
    directory_name = os.path.dirname(file_name)
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)
    output_file = open(file_name, "wt")
    op_gen_lines = []
    for operand in operands:
        replacement = {"name": operand.name}
        if operand.has_property(properties.SYMMETRIC):
            replacement["shape"] = "Symmetric"
        elif operand.has_property(properties.DIAGONAL):
            replacement["shape"] = "Diagonal"
        elif operand.has_property(properties.LOWER_TRIANGULAR):
            replacement["shape"] = "Triangular"
        else:
            replacement["shape"] = "MISSING"
        replacement["size"] = ", ".join(str(val) for val in operand.size)

        replacement["properties"] = ""

        op_gen_lines.append(op_gen_line_template.format(**replacement))

        # print(operand)
    op_gen_lines = "\n".join(op_gen_lines)
    # print(op_gen_lines)
    op_gen_file = op_gen_file_template.format(textwrap.indent(op_gen_lines, "    "), output_str)
    output_file.write(op_gen_file)
    output_file.close()
