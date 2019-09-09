from ... import config

from . import operand_generation, runner, reference_code


def generate_experiment_code(name, equations, algorithm_name, num_threads, k_best=False, k=1):

    reference_code.generate_reference_code(name, equations)

    operand_generation.generate_operand_generator(name, equations)

    algorithms = [("experiments", algorithm_name.format(0))]
    runner.generate_runner(name, algorithms, num_threads)

    if k_best:
        k_best_algorithms = [("k_best", algorithm_name.format(i)) for i in range(k)]
        for threads in num_threads:
            runner.runner_to_file("runner_k_best", name, config.Language.Julia, threads, k_best_algorithms)
        