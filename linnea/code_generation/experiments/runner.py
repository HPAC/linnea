
import textwrap
import os.path

from ... import config
from .. import utils

algorithm_inclusion = {config.Language.Julia: """include("{0}/{1}.jl")""",
                    config.Language.Cpp: ""}

algorithm_tmp_result = {config.Language.Julia: """result_{0} = collect({0}(map(MatrixGenerator.unwrap, matrices)...))""",
                    config.Language.Cpp: ""}

algorithm_test = {config.Language.Julia: "@test isapprox(result_{0}, result_recommended)",
                    config.Language.Cpp: ""}

algorithm_plot = {config.Language.Julia: """Benchmarker.add_data(plotter, ["{0}"], Benchmarker.measure(20, {0}, map(MatrixGenerator.unwrap, matrices)...) );""",
                    config.Language.Cpp: ""}



def runner_to_file(output_name, language, algorithms=[]):

    if language is config.Language.Julia:
        file_name = "runner.jl"
    elif language is config.Language.Cpp:
        file_name = "runner.cpp"
    elif language is config.Language.Matlab:
        file_name = "runner.m"
    else:
        raise config.LanguageOptionNotImplemented()

    file_path = os.path.join(config.output_path, output_name, language.name, file_name)
    output_file = open(file_path, "wt", encoding='utf-8')
    inclusions = []
    tmp_results = []
    tests = []
    plots = []
    incl_format = algorithm_inclusion.get(language)
    tmp_format = algorithm_tmp_result.get(language)
    test_format = algorithm_test.get(language)
    plot_format = algorithm_plot.get(language)

    for subdir_name, algorithm_name in algorithms:
        inclusions.append(incl_format.format(subdir_name, algorithm_name))
        tmp_results.append(tmp_format.format(algorithm_name))
        tests.append(test_format.format(algorithm_name))
        plots.append(plot_format.format(algorithm_name))


    runner_template = utils.get_template(file_name, language)

    if language == config.Language.Julia:
        runner_file = runner_template.format(
            "\n".join(inclusions),
            "\n".join(tmp_results),
            "\n".join(tests),
            output_name,
            "\n".join(plots),
        )
    else:
        runner_file = runner_template.format(output_name)

    output_file.write(runner_file)
    output_file.close()


def generate_cmake_script(output_name):

    file_name = "CMakeLists.txt"
    script_template = utils.get_template(file_name, config.Language.Cpp)
    file_path = os.path.join(config.output_path, output_name, config.Language.Cpp.name, file_name)
    output_file = open(file_path, "wt")
    output_file.write(script_template.format(output_name))
    output_file.close()


def generate_runner(output_name, algorithms):

    runner_to_file(output_name, algorithms=algorithms, language=config.Language.Julia)
    runner_to_file(output_name, language=config.Language.Matlab)
    runner_to_file(output_name, language=config.Language.Cpp)
    generate_cmake_script(output_name)
