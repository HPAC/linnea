
import os.path
import textwrap

from ... import config
from .. import utils

from ...algebra.properties import Property


operands_mapping_julia = {Property.SYMMETRIC: "Shape.Symmetric",
                          Property.DIAGONAL: "Shape.Diagonal",
                          Property.LOWER_TRIANGULAR: "Shape.LowerTriangular",
                          Property.UPPER_TRIANGULAR: "Shape.UpperTriangular",
                          Property.SPD: "Properties.SPD",
                          Property.SPSD: "Properties.SPD", # TODO does this make sense?
                          Property.ORTHOGONAL: "Properties.Orthogonal"
                         }

operands_mapping_matlab = {Property.SYMMETRIC: "Shape.Symmetric()",
                          Property.DIAGONAL: "Shape.Diagonal()",
                          Property.LOWER_TRIANGULAR: "Shape.LowerTriangular()",
                          Property.UPPER_TRIANGULAR: "Shape.UpperTriangular()",
                          Property.SPD: "Properties.SPD()",
                          Property.SPSD: "Properties.SPD()", # TODO does this make sense?
                          Property.ORTHOGONAL: "Properties.Orthogonal()"
                         }

operands_mapping_cpp = {Property.SYMMETRIC: "generator::shape::self_adjoint{}",
                        Property.DIAGONAL: "generator::shape::diagonal{}",
                        Property.LOWER_TRIANGULAR: "generator::shape::lower_triangular{}",
                        Property.UPPER_TRIANGULAR: "generator::shape::upper_triangular{}",
                        Property.SPD: "generator::property::spd{}",
                        Property.SPSD: "generator::property::spd{}", # TODO does this make sense?
                        Property.ORTHOGONAL: "generator::property::orthogonal{}"
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
        op_gen_line_template = "{name}::{type} = generate(({size}), [{properties}])"
        # random_operand = ["Properties.Random(20, 21)"]
        def random_operand(l, u):
            return "Properties.Random({}, {})".format(l, u)

    elif language is config.Language.Cpp:
        file_name = "operand_generator.hpp"
        op_gen_line_template = "auto {name} = gen.generate({{{size}}}, {properties});"
        # random_operand = ["generator::property::random{}"]
        def random_operand(l, u):
            return "generator::property::random{}"

    elif language is config.Language.Matlab:
        file_name = "operand_generator.m"
        op_gen_line_template = "{name} = generate([{size}], {properties});"
        # random_operand = ["Properties.Random([20, 21])"]
        def random_operand(l, u):
            return "Properties.Random([{}, {}])".format(l, u)


    else:
        raise config.LanguageOptionNotImplemented()


    file_path = os.path.join(config.output_code_path, output_name, language.name, file_name)
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

        replacement["type"] = utils.operand_type(operand, True)

        # Special case - scalar generation for C++
        if language == config.Language.Cpp and operand.has_property(Property.SCALAR):
            op_gen_lines.append("{0} {1} = std::uniform_real_distribution<{0}>()(gen.rng());".format(
                "double" if config.float64 else "float",
                operand.name
            ))
            continue

        property_replacements = []
        for prop in [Property.SYMMETRIC, Property.DIAGONAL, Property.LOWER_TRIANGULAR,
                     Property.UPPER_TRIANGULAR, Property.ORTHOGONAL]:
            if prop in operand.properties:
                property_replacements.append(map_operand(language, prop))
        if Property.SPD in operand.properties or Property.SPSD in operand.properties:
            property_replacements.append(map_operand(language, Property.SPD))
        if not any((prop in operand.properties) for prop in [Property.SYMMETRIC, Property.DIAGONAL, Property.LOWER_TRIANGULAR, Property.UPPER_TRIANGULAR]):
            if language == config.Language.Julia:
                property_replacements.append("Shape.General")
            elif language == config.Language.Matlab:
                property_replacements.append("Shape.General()")

        if not Property.SPD in operand.properties and not Property.SPSD in operand.properties and not Property.ORTHOGONAL in operand.properties:
            if Property.DIAGONAL in operand.properties:
                property_replacements.append(random_operand(10, 11))
            elif Property.LOWER_TRIANGULAR in operand.properties:
                property_replacements.append(random_operand(10, 11))
            elif Property.UPPER_TRIANGULAR in operand.properties:
                property_replacements.append(random_operand(10, 11))
            elif Property.SYMMETRIC in operand.properties:
                property_replacements.append(random_operand(10, 11))
            elif operand.has_property(Property.SCALAR):
                property_replacements.append(random_operand(0.5, 1.5))
            else:
                property_replacements.append(random_operand(-1, 1))

        if language == config.Language.Cpp:
            if Property.VECTOR in operand.properties:
                if operand.size[0] == 1:
                    property_replacements.append("generator::shape::row_vector{}")
                else:
                    property_replacements.append("generator::shape::col_vector{}")
            if operand.size[0] != operand.size[1]:
                property_replacements.append("generator::shape::not_square{}")
                
        replacement["size"] = ",".join(map(str, operand.size))

        replacement["properties"] = ", ".join(property_replacements)
        op_gen_lines.append(op_gen_line_template.format(**replacement))

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

