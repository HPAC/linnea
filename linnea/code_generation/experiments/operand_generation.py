
import os.path
import textwrap

from ... import config
from .. import utils

from ...algebra.properties import Property as properties


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


def map_operand(language, property):
    return operands_mapping[language][property]


# TODO: somehow refactor this. it's ugly
def operand_generator_to_file(output_name, operands, output_str, language = config.language):

    if language is config.Language.Julia:
        file_name = "operand_generator.jl"
        op_gen_line_template = "{name} = generate(({size}), [{properties}])"
        random_operand = "Properties.Random", "Shape.General"

    elif language is config.Language.Cpp:
        file_name = "operand_generator.hpp"
        op_gen_line_template = "auto {name} = gen.generate({{{size}}}, {properties});"
        random_operand = "generator::property::random{}"

    elif language is config.Language.Matlab:
        file_name = "operand_generator.m"
        op_gen_line_template = "{name} = generate([{size}], {properties});"
        random_operand = "Properties.Random()", "Shape.General()"

    else:
        raise config.LanguageOptionNotImplemented()


    file_path = os.path.join(config.output_path, output_name, language.name, file_name)
    directory_name = os.path.dirname(file_path)
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)
    output_file = open(file_path, "wt")
    op_gen_lines = []

    for idx, operand in enumerate(operands, 1):
        if language != config.Language.Matlab:
            replacement = {"name": operand.name}
        else:
            replacement = {"name": "out{{ {0} }}".format(idx)}

        # Special case - scalar generation for C++
        if language == config.Language.Cpp and properties.SCALAR in operand.properties:
            op_gen_lines.append("{0} {1} = std::uniform_real_distribution<{0}>()(gen.rng());".format(
                "double" if config.float64 else "float",
                operand.name
            ))
            continue

        property_replacements = []
        for prop in [properties.SYMMETRIC, properties.DIAGONAL, properties.LOWER_TRIANGULAR,
                     properties.UPPER_TRIANGULAR, properties.SPD]:
            if prop in operand.properties:
                property_replacements.append(map_operand(language, prop))
        if not properties.SPD in operand.properties:
            property_replacements.extend(random_operand)
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
        op_gen_lines.append(op_gen_line_template.format(**replacement))

        # print(operand)
    op_gen_lines = "\n".join(op_gen_lines)
    op_gen_file_template = utils.get_template(file_name, language)
    op_gen_file = op_gen_file_template.format(textwrap.indent(op_gen_lines, "    "), output_str)
    output_file.write(op_gen_file)
    output_file.close()
    if config.verbosity >= 2:
        print("Generate operand generator file {}".format(file_path))

def generate_operand_generator(output_name, equations):
    input, output = equations.input_output()
    input_str = ", ".join([operand.name for operand in input])

    operand_generator_to_file(output_name, input, input_str)
    operand_generator_to_file(output_name, input, input_str, language=config.Language.Cpp)
    operand_generator_to_file(output_name, input, input_str, language=config.Language.Matlab)

