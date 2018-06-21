
import os.path
import textwrap

from ... import config

from ...algebra.properties import Property as properties

filename_extension = {config.Language.Julia: ".jl",
                      config.Language.Cpp: ".hpp",
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

op_gen_line_template = {config.Language.Julia : "{name} = generate(({size}), [{properties}])",
                        config.Language.Cpp: "auto {name} = gen.generate({{{size}}}, {properties});",
                        config.Language.Matlab: "{name} = generate([{size}], {properties});"}

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


def map_operand(language, property):
    return operands_mapping.get(language).get(property)


# TODO: somehow refactor this. it's ugly
def operand_generator_to_file(output_name, operands, output_str, language = config.language, name_addition = ""):
    additional_name = "_{}".format(name_addition) if name_addition else ""
    file_name = os.path.join(config.output_path, output_name, language.name,
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

        # Special case - scalar generation for C++
        if language == config.Language.Cpp and properties.SCALAR in operand.properties:
            op_gen_lines.append("{0} {1} = std::uniform_real_distribution<{0}>()(gen.rng());".format(
                "double" if config.float64 else "float",
                operand.name
            ))
            counter += 1
            continue

        property_replacements = []
        for prop in [properties.SYMMETRIC, properties.DIAGONAL, properties.LOWER_TRIANGULAR,
                     properties.UPPER_TRIANGULAR, properties.SPD]:
            if prop in operand.properties:
                property_replacements.append(map_operand(language, prop))
        if not properties.SPD in operand.properties:
            property_replacements.extend(random_operand.get(language))
        if language == config.Language.Cpp:
            if properties.VECTOR in operand.properties:
                if operand.size[0] == 1:
                    property_replacements.append("generator::shape::row_vector{}")
                else:
                    property_replacements.append("generator::shape::col_vector{}")
            if operand.size[0] != operand.size[1]:
                property_replacements.append("generator::shape::not_square{}")
                
        replacement["size"] = ",".join(map(str, operand.size))

        replacement["properties"] = ", ".join(property_replacements)
        op_gen_lines.append(op_gen_line_template.get(language).format(**replacement))

        # print(operand)
        counter += 1
    op_gen_lines = "\n".join(op_gen_lines)
    op_gen_file = op_gen_file_template.get(language).format(textwrap.indent(op_gen_lines, "    "), output_str)
    output_file.write(op_gen_file)
    output_file.close()
    if config.verbosity >= 2:
        print("Generate operand generator file {}".format(file_name))

def generate_operand_generator(output_name, equations):
    input, output = equations.input_output()
    input_str = ", ".join([operand.name for operand in input])

    operand_generator_to_file(output_name, input, input_str)
    operand_generator_to_file(output_name, input, input_str, language=config.Language.Cpp)
    operand_generator_to_file(output_name, input, input_str, language=config.Language.Matlab)

