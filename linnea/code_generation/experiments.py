
from ..algebra.properties import Property as properties
from .. import config
import textwrap
import os.path


op_gen_file_template = textwrap.dedent(
                        """
                        push!(LOAD_PATH, "/Users/henrik/Promotion/code/linalg_tests/julia/")
                        using Generator

                        function operand_generator()
                        {}
                            return ({})
                        end
                        """)

op_gen_line_template = """{name} = generate({size}, [{properties}])"""

experiment_template = textwrap.dedent(
                            """
                            """
                            )


def operand_generator_to_file(output_name, operands, output_str):
    file_name = os.path.join(config.output_path, config.language.name, output_name, "operand_generator.jl")
    directory_name = os.path.dirname(file_name)
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)
    output_file = open(file_name, "wt")
    op_gen_lines = []
    for operand in operands:
        replacement = {"name": operand.name}
        property_replacements = []
        if properties.SYMMETRIC in operand.properties:
            property_replacements.append("Shape.Symmetric")
        if properties.DIAGONAL in operand.properties:
            property_replacements.append("Shape.Diagonal")
        if properties.LOWER_TRIANGULAR in operand.properties:
            property_replacements.append("Shape.LowerTriangular")
        if properties.UPPER_TRIANGULAR in operand.properties:
            property_replacements.append("Shape.UpperTriangular")
        if properties.SPD in operand.properties:
            property_replacements.append("Properties.SPD")
        else:
            property_replacements.append("Properties.Random")
            property_replacements.append("Shape.General")


        replacement["size"] = str(operand.size)

        replacement["properties"] = ", ".join(property_replacements)

        op_gen_lines.append(op_gen_line_template.format(**replacement))

        # print(operand)
    op_gen_lines = "\n".join(op_gen_lines)
    # print(op_gen_lines)
    op_gen_file = op_gen_file_template.format(textwrap.indent(op_gen_lines, "    "), output_str)
    output_file.write(op_gen_file)
    output_file.close()
