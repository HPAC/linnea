
import textwrap
import os.path

from ... import config
from .. import utils

algorithm_inclusion = {config.Language.Julia: """include("{0}/{1}.jl")""",
                         config.Language.Cpp: ""}

algorithm_test = {config.Language.Julia: textwrap.dedent(
                                         """\
                                         result = collect({0}(map(MatrixGenerator.unwrap, map(copy, matrices))...)[1])
                                         test_result = isapprox(result, result_recommended, rtol=1e-3)
                                         @test test_result # this somehow avoids too much memory being used\
                                         """),
                    config.Language.Cpp: ""}

algorithm_plot = {config.Language.Julia: """Benchmarker.add_data(plotter, ["{0}"; {1}], Benchmarker.measure({2}, {0}, map(MatrixGenerator.unwrap, matrices)...) );""",
                    config.Language.Cpp: ""}



def runner_to_file(runner_name, output_name, language, num_threads, repetitions, algorithms=[]):
    
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
    tests = []
    plots = []
    incl_format = algorithm_inclusion.get(language)
    test_format = algorithm_test.get(language)
    plot_format = algorithm_plot.get(language)

    for subdir_name, algorithm_name in algorithms:
        includes.append(incl_format.format(subdir_name, algorithm_name))
        tests.append(test_format.format(algorithm_name))
        plots.append(plot_format.format(algorithm_name, num_threads, repetitions))


    runner_template = utils.get_template(template_name, language)

    if language == config.Language.Julia:
        runner_file = runner_template.format(
            num_threads = num_threads,
            repetitions = repetitions,
            includes = "\n".join(includes),
            tests = textwrap.indent("\n".join(tests), "    "),
            experiment_name = output_name,
            measurements = textwrap.indent("\n".join(plots), "    "),
        )
    else:
        runner_file = runner_template.format(num_threads = num_threads, repetitions = repetitions, experiment_name = output_name)

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


def generate_runner(output_name, algorithms, num_threads, repetitions):

    for threads in num_threads:
        runner_to_file("runner", output_name, config.Language.Julia, threads, repetitions, algorithms=algorithms)
        runner_to_file("runner", output_name, config.Language.Matlab, threads, repetitions)
        runner_to_file("runner", output_name, config.Language.Cpp, threads, repetitions)
    generate_cmake_script(output_name, num_threads)
