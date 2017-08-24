
from ..algebra.properties import Property as properties
from .. import config
import textwrap
import os.path


op_gen_file_template = {config.Language.Julia: textwrap.dedent(
                            """
                            using Generator

                            function operand_generator()
                            {}
                                return ({})
                            end
                            """
                        ),
                        config.Language.Cpp: textwrap.dedent(
                            """
                            #include <generator/generator.hpp>

                            template<typename T, typename Gen>
                            decltype(auto) operand_generator()
                            {{
                            {}
                                return std::make_tuple({});
                            }}
                            """
                        )
                        }
op_gen_line_template = {config.Language.Julia : "{name} = generate({size}, [{properties}])",
                        config.Language.Cpp: "auto {name} = generator::generate<T, Gen>({size}, {properties});"}

experiment_template = textwrap.dedent(
                            """
                            """
                            )

filename_extension = {config.Language.Julia : ".jl",
                      config.Language.Cpp: ".hpp"}


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

        op_gen_lines.append(op_gen_line_template.get(language).format(**replacement))

        # print(operand)
    op_gen_lines = "\n".join(op_gen_lines)
    # print(op_gen_lines)
    op_gen_file = op_gen_file_template.get(language).format(textwrap.indent(op_gen_lines, "    "), output_str)
    output_file.write(op_gen_file)
    output_file.close()
