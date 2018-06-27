from .. import utils
from ... import config

import os.path

def generate_reference_code(output_name, equations):
    input, output = equations.input_output()
    input_str = ", ".join([operand.name for operand in input])
    output_str = ", ".join([operand.name for operand in output])

    utils.remove_files(os.path.join(config.output_path, output_name, config.Language.Julia.name, "reference"))
    utils.remove_files(os.path.join(config.output_path, output_name, config.Language.Cpp.name, "reference"))
    utils.remove_files(os.path.join(config.output_path, output_name, config.Language.Matlab.name, "reference"))

    utils.algorithm_to_file(output_name, "reference", "naive_blaze",
                        equations.to_cpp_expression(config.CppLibrary.Blaze),
                        input_str, output_str,
                        config.Language.Cpp, ".hpp")
    utils.algorithm_to_file(output_name, "reference", "naive_eigen",
                        equations.to_cpp_expression(config.CppLibrary.Eigen),
                        input_str, output_str,
                        config.Language.Cpp, ".hpp")
    utils.algorithm_to_file(output_name, "reference", "naive_armadillo",
                        equations.to_cpp_expression(config.CppLibrary.Armadillo),
                        input_str, output_str,
                        config.Language.Cpp, ".hpp")
    utils.algorithm_to_file(output_name, "reference", "naive",
                        equations.to_julia_expression(),
                        input_str, output_str,
                        config.Language.Julia)
    utils.algorithm_to_file(output_name, "reference", "naive",
                        equations.to_matlab_expression(), input_str, output_str,
                        config.Language.Matlab, ".m")


    equations = utils.replace_linsolve_left(equations)

    # those libraries do not support solving linear systems A*inv(B)
    utils.algorithm_to_file(output_name, "reference", "recommended_eigen",
                        equations.to_cpp_expression(config.CppLibrary.Eigen, recommended=True),
                        input_str, output_str,
                        config.Language.Cpp, ".hpp")
    utils.algorithm_to_file(output_name, "reference", "recommended_armadillo",
                        equations.to_cpp_expression(config.CppLibrary.Armadillo, recommended=True),
                        input_str, output_str,
                        config.Language.Cpp, ".hpp")
    
    equations = utils.replace_linsolve_right(equations)

    utils.algorithm_to_file(output_name, "reference", "recommended",
                        equations.to_julia_expression(recommended=True),
                        input_str, output_str,
                        config.Language.Julia)
    utils.algorithm_to_file(output_name, "reference", "recommended",
                        equations.to_matlab_expression(recommended=True),
                        input_str, output_str,
                        config.Language.Matlab, ".m")