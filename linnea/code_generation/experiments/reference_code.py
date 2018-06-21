from .. import utils

from ... import config


def generate_reference_code(output_name, equations):
    input, output = equations.input_output()
    input_str = ", ".join([operand.name for operand in input])
    output_str = ", ".join([operand.name for operand in output])

    utils.algorithm_to_file(output_name, "reference", "naive", equations.to_julia_expression(), input_str, output_str, config.Language.Julia)
    utils.algorithm_to_file(output_name, "reference", "recommended", equations.to_julia_expression(recommended=True),
                          input_str, output_str, config.Language.Julia)
    utils.algorithm_to_file(output_name, "reference", "naive", equations.to_cpp_expression(config.CppLibrary.Blaze),
                          input_str, output_str, config.Language.Cpp, ".hpp", "blaze")
    utils.algorithm_to_file(output_name, "reference", "naive", equations.to_cpp_expression(config.CppLibrary.Eigen),
                          input_str, output_str, config.Language.Cpp, ".hpp", "eigen")
    utils.algorithm_to_file(output_name, "reference", "naive", equations.to_cpp_expression(config.CppLibrary.Armadillo),
                          input_str, output_str, config.Language.Cpp, ".hpp", "armadillo")
    utils.algorithm_to_file(output_name, "reference", "recommended",
                          equations.to_cpp_expression(config.CppLibrary.Eigen, recommended=True),
                          input_str, output_str, config.Language.Cpp, ".hpp", "eigen")
    utils.algorithm_to_file(output_name, "reference", "recommended",
                          equations.to_cpp_expression(config.CppLibrary.Armadillo, recommended=True),
                          input_str, output_str, config.Language.Cpp, ".hpp", "armadillo")
    utils.algorithm_to_file(output_name, "reference", "naive", equations.to_matlab_expression(), input_str, output_str,
                          config.Language.Matlab, ".m")
    utils.algorithm_to_file(output_name, "reference", "recommended", equations.to_matlab_expression(recommended=True),
                          input_str, output_str, config.Language.Matlab, ".m")