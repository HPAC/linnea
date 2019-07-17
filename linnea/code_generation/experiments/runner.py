
import textwrap
import os.path

from ... import config
from .. import utils

algorithm_inclusion = {config.Language.Julia: """include("{0}/{1}.jl")""",
                    config.Language.Cpp: ""}

algorithm_tmp_result = {config.Language.Julia: """result_{0} = collect({0}(map(MatrixGenerator.unwrap, map(copy, matrices))...)[1])""",
                    config.Language.Cpp: ""}

algorithm_test = {config.Language.Julia: "@test isapprox(result_{0}, result_recommended, rtol=1e-3)",
                    config.Language.Cpp: ""}

algorithm_plot = {config.Language.Julia: """Benchmarker.add_data(plotter, ["{0}"; {1}], Benchmarker.measure(20, {0}, map(MatrixGenerator.unwrap, matrices)...) );""",
                    config.Language.Cpp: ""}



def runner_to_file(runner_name, output_name, language, num_threads, algorithms=[]):

    if language is config.Language.Julia:
        file_name = "{}_t{}.jl".format(runner_name, num_threads)
        template_name = "runner.jl"
    elif language is config.Language.Cpp:
        file_name = "{}_t{}.cpp".format(runner_name, num_threads)
        template_name = "runner.cpp"
    elif language is config.Language.Matlab:
        file_name = "{}_t{}.m".format(runner_name, num_threads)
        template_name = "runner.m"
    else:
        raise config.LanguageOptionNotImplemented()

    file_path = os.path.join(config.output_code_path, output_name, language.name, file_name)
    output_file = open(file_path, "wt", encoding='utf-8')
    includes = []
    tmp_results = []
    tests = []
    plots = []
    incl_format = algorithm_inclusion.get(language)
    tmp_format = algorithm_tmp_result.get(language)
    test_format = algorithm_test.get(language)
    plot_format = algorithm_plot.get(language)

    for subdir_name, algorithm_name in algorithms:
        includes.append(incl_format.format(subdir_name, algorithm_name))
        tmp_results.append(tmp_format.format(algorithm_name))
        tests.append(test_format.format(algorithm_name))
        plots.append(plot_format.format(algorithm_name, num_threads))


    runner_template = utils.get_template(template_name, language)

    if language == config.Language.Julia:
        runner_file = runner_template.format(
            num_threads = num_threads,
            includes = "\n".join(includes),
            tmp_results = textwrap.indent("\n".join(tmp_results), "    "),
            tests = textwrap.indent("\n".join(tests), "    "),
            experiment_name = output_name,
            measurements = textwrap.indent("\n".join(plots), "    "),
        )
    else:
        runner_file = runner_template.format(num_threads = num_threads, experiment_name = output_name)

    if config.verbosity >= 2:
        print("Generate runner file {}".format(file_path))
    output_file.write(runner_file)
    output_file.close()


def generate_cmake_script(output_name, num_threads):

    file_name = "CMakeLists.txt"
    script_template = utils.get_template(file_name, config.Language.Cpp)
    executables = []
    for threads in num_threads:
        executables.append("add_executable(runner_t{0} runner_t{0}.cpp)\ntarget_link_libraries(runner_t{0} PRIVATE libtests)\nset_target_properties(runner_t{0} PROPERTIES CXX_STANDARD 14)\n".format(threads))
    script = script_template.format(experiment_name = output_name, executables = "\n".join(executables))
    file_path = os.path.join(config.output_code_path, output_name, config.Language.Cpp.name, file_name)
    output_file = open(file_path, "wt")
    output_file.write(script)
    output_file.close()


def generate_runner(output_name, algorithms, num_threads):

    for threads in num_threads:
        runner_to_file("runner", output_name, config.Language.Julia, threads, algorithms=algorithms)
        runner_to_file("runner", output_name, config.Language.Matlab, threads)
        runner_to_file("runner", output_name, config.Language.Cpp, threads)
    generate_cmake_script(output_name, num_threads)
