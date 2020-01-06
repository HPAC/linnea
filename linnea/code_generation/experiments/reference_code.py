from .. import utils
from ... import config

import os.path
import textwrap

def equations_to_code(equations, input, output, function_name, language, library=None):
    if language == config.Language.Julia:
        template = utils.get_template("algorithm_experiments.jl", language)
        input_str = ", ".join(["{}::{}".format(operand.name, utils.operand_type(operand, True)) for operand in input])
        output_str = ", ".join([operand.name for operand in output])
        code = equations.to_julia_expression()
        algorithm_str = template.format(function_name, input_str, textwrap.indent(code, "    "), output_str)
    elif language == config.Language.Matlab:
        template = utils.get_template("algorithm.m", language)
        input_str = ", ".join([operand.name for operand in input])
        output_str = ", ".join([operand.name for operand in output])
        code = equations.to_matlab_expression() 
        algorithm_str = template.format(function_name, input_str, textwrap.indent(code, "    "), output_str)
    elif language == config.Language.Cpp:
        template = utils.get_template("algorithm.cpp", language)
        types_str = ", ".join("typename Type_{}".format(op) for op in input)
        args_str = ", ".join("Type_{0} && {0}".format(op) for op in input)
        output_len = len(output)
        output_str = ", ".join([operand.name for operand in output])
        if output_len > 1:
            output_string = "return std::make_tuple({});".format(output_str)
        else:
            ret_string = textwrap.dedent(
                         """\
                         typedef std::remove_reference_t<decltype({})> return_t;
                         return return_t({});\
                         """)
            output_string = ret_string.format(output_str, output_str)
        code = equations.to_cpp_expression(library)
        algorithm_str = template.format(function_name, types_str, args_str, textwrap.indent(code, "    "), textwrap.indent(output_string, "    "))
    else:
        raise config.LanguageOptionNotImplemented()
    return algorithm_str


def generate_reference_code(output_name, equations):
    input, output = equations.input_output()
    input_str = ", ".join([operand.name for operand in input])
    output_str = ", ".join([operand.name for operand in output])

    julia_input_str = ", ".join(["{}::{}".format(operand.name, utils.operand_type(operand, True)) for operand in input])

    output_path = os.path.join(config.output_code_path, output_name)

    julia_path = os.path.join(output_path, config.Language.Julia.name, "reference")
    cpp_path = os.path.join(output_path, config.Language.Cpp.name, "reference")
    matlab_path = os.path.join(output_path, config.Language.Matlab.name, "reference")

    utils.remove_files(julia_path)
    utils.remove_files(cpp_path)
    utils.remove_files(matlab_path)

    file_name = os.path.join(cpp_path, "naive_eigen.hpp")
    code = equations_to_code(equations, input, output, "naive_eigen", config.Language.Cpp, config.CppLibrary.Eigen)
    utils.to_file(file_name, code)

    file_name = os.path.join(cpp_path, "naive_armadillo.hpp")
    code = equations_to_code(equations, input, output, "naive_armadillo", config.Language.Cpp, config.CppLibrary.Armadillo)
    utils.to_file(file_name, code)

    file_name = os.path.join(julia_path, "naive.jl")
    code = equations_to_code(equations, input, output, "naive", config.Language.Julia)
    utils.to_file(file_name, code)

    file_name = os.path.join(matlab_path, "naive.m")
    code = equations_to_code(equations, input, output, "naive", config.Language.Matlab)
    utils.to_file(file_name, code)


    equations = utils.replace_linsolve_left(equations)

    # those libraries do not support solving linear systems A*inv(B)
    file_name = os.path.join(cpp_path, "recommended_eigen.hpp")
    code = equations_to_code(equations, input, output, "recommended_eigen", config.Language.Cpp, config.CppLibrary.Eigen)
    utils.to_file(file_name, code)

    file_name = os.path.join(cpp_path, "recommended_armadillo.hpp")
    code = equations_to_code(equations, input, output, "recommended_armadillo", config.Language.Cpp, config.CppLibrary.Armadillo)
    utils.to_file(file_name, code)

    
    equations = utils.replace_linsolve_right(equations)

    file_name = os.path.join(julia_path, "recommended.jl")
    code = equations_to_code(equations, input, output, "recommended", config.Language.Julia)
    utils.to_file(file_name, code)

    file_name = os.path.join(matlab_path, "recommended.m")
    code = equations_to_code(equations, input, output, "recommended", config.Language.Matlab)
    utils.to_file(file_name, code)
    