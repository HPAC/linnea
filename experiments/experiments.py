import pandas as pd
import time
import numpy
import itertools
import argparse
import collections
import random
import os
import logging

logging.basicConfig(level=logging.INFO,
                    format='%(name)-2s: %(levelname)-2s %(message)s')

import linnea.config

linnea.config.init()

from linnea.config import Strategy, Language
import linnea.examples.examples
import linnea.examples.PLDI as PLDI
from linnea.derivation.graph.constructive import DerivationGraph as CDGraph
from linnea.derivation.graph.exhaustive import DerivationGraph as EDGraph
from linnea.code_generation.experiments import operand_generation, runner, reference_code

from random_expressions import generate_equation
from jobscripts import generate_scripts

def measure(example, name, strategy, merge, reps=10):

    times = []

    if strategy is Strategy.constructive:
        DerivationGraph = CDGraph
    elif strategy is Strategy.exhaustive:
        DerivationGraph = EDGraph
    else:
        raise NotImplementedError()

    for i in range(reps):
        linnea.config.clear_all()
        if hasattr(example, "init"):
            # calls initialization that have to be done before each repetition
            example.init()
        graph = DerivationGraph(example.eqns)
        t_start = time.perf_counter()
        trace = graph.derivation(
                            solution_nodes_limit=100,
                            iteration_limit=15,
                            merging=merge,
                            dead_ends=True)
        graph.write_output(code=True,
                           derivation=False,
                           output_name=name,
                           experiment_code=False,
                           algorithms_limit=1,
                           graph=False)
        t_end = time.perf_counter()
        times.append(t_end-t_start)
    data = [numpy.mean(times),
            numpy.std(times),
            numpy.min(times),
            numpy.max(times),
            graph.nodes_count(),
            bool(graph.terminal_nodes())]
    return data

def generate(experiment, example, name, strategy):

    strategy_str = None 
    if strategy is Strategy.constructive:
        DerivationGraph = CDGraph
        algorithm_name = "algorithm{}c"
        strategy_str = "c"
    elif strategy is Strategy.exhaustive:
        DerivationGraph = EDGraph
        algorithm_name = "algorithm{}e"
        strategy_str = "e"
    else:
        raise NotImplementedError()

    linnea.config.clear_all()
    if hasattr(example, "init"):
        # calls initialization that have to be done before each repetition
        example.init()

    print(example.eqns)

    graph = DerivationGraph(example.eqns)
    trace = graph.derivation(
                        solution_nodes_limit=100,
                        iteration_limit=15,
                        merging=True,
                        dead_ends=True)
    graph.write_output(code=True,
                       derivation=False,
                       output_name=name,
                       experiment_code=False,
                       algorithms_limit=1,
                       graph=False,
                       subdir_name=strategy.name,
                       algorithm_name=algorithm_name)

    _, cost = graph.shortest_path()
    data = example.eqns.get_data()
    file_path = os.path.join(linnea.config.results_path, experiment, "intensity", strategy_str, name + "_intensity" + ".csv")
    mindex = pd.MultiIndex.from_arrays([[name], [algorithm_name.format(0)]], names=("example", "algorithm"))
    dframe = pd.DataFrame([[data, cost, cost/data]], index=mindex, columns=["data", "cost", "intensity"])
    dframe.to_csv(file_path)

def main():

    parser = argparse.ArgumentParser(prog="experiments")
    parser.add_argument("mode", choices=["time_generation", "generate_code", "jobscripts"])
    parser.add_argument("experiment", choices=["random", "pldi"])
    parser.add_argument("-m", "--merging", choices=["true", "false", "both"], default="true")
    parser.add_argument("-j", "--jobindex", help="Job index.", type=int, default=0)
    parser.add_argument("-r", "--repetitions", help="Number of repetitions.", type=int)
    parser.add_argument("-c", "--constructive", action="store_true", help="Use constructive strategy.")
    parser.add_argument("-e", "--exhaustive", action="store_true", help="Use exhaustive strategy.")
    parser.add_argument("-f", "--reference", action="store_true", help="Generate reference code.")
    parser.add_argument("-l", "--config", type=str, default=None, help="Specify configuration file.")
    args = parser.parse_args()

    PLDI_examples = [
        PLDI.Example01(),
        PLDI.Example02(),
        PLDI.Example03(),
        PLDI.Example04(),
        PLDI.Example05(),
        PLDI.Example06(),
        PLDI.Example07(),
        PLDI.Example08(),
        PLDI.Example09(),
        PLDI.Example10(),
        PLDI.Example11(),
        PLDI.Example12(),
        PLDI.Example13(),
        PLDI.Example14(),
        PLDI.Example15(),
        PLDI.Example16(),
        PLDI.Example17(),
        PLDI.Example18(),
        PLDI.Example19(),
        PLDI.Example20(),
        PLDI.Example21(),
        PLDI.Example22(),
        PLDI.Example23(),
        PLDI.Example24(),
        PLDI.Example25(),
    ]

    ExampleContainer = collections.namedtuple("ExampleContainer", ["eqns"])

    random.seed(0)
    rand_exprs = [ExampleContainer(generate_equation(random.randint(4, 7))) for _ in range(100)]

    # also add init for new examples

    # TODO write parameters (repetitions) to file
    # TODO write description of experiment (name, operand sizes, equations) to file

    if args.config:
        linnea.config.load_config(config_file=args.config)

    if args.experiment == "random":
        examples = rand_exprs
    elif args.experiment == "pldi":
        examples = PLDI_examples
    else:
        return


    JobExample = collections.namedtuple("JobExample", ["example", "name"])

    if args.jobindex == 0:
        job_examples = []
        for idx, example in enumerate(examples, 1):
            name = "{}{:03}".format(args.experiment, idx)
            job_examples.append(JobExample(example, name))
    else:
        name = "{}{:03}".format(args.experiment, args.jobindex)
        job_examples = [JobExample(examples[args.jobindex-1], name)]


    strategies = []
    algorithms = []
    strategy_strs = []
    if args.constructive:
        strategies.append(Strategy.constructive)
        algorithms.append(("constructive", "algorithm0c"))
        strategy_strs.append("c")
    if args.exhaustive:
        strategies.append(Strategy.exhaustive)
        algorithms.append(("exhaustive", "algorithm0e"))
        strategy_strs.append("e")

    if args.mode == "time_generation":

        if not strategies:
            return

        linnea.config.set_verbosity(0)

        merging_args = []
        merging_labels = []
        if args.merging == "true" or args.merging == "both":
            merging_args.append(True)
            merging_labels.append("merging")
        elif args.merging == "false" or args.merging == "both":
            merging_args.append(False)
            merging_labels.append("no_merging")

        mindex = pd.MultiIndex.from_product([
            [example.name for example in job_examples],
            [strategy.name for strategy in strategies],
            merging_labels],
            names=["example", "strategy", "merging"])

        col_index = pd.Index(["mean", "std", "min", "max", "nodes", "solution"])

        if args.jobindex == 0:
            output_file = "generation.csv"
        else:
            output_file = "generation{:03}.csv".format(args.jobindex)

        data = []
        for example, name in job_examples:
            for strategy, merge in itertools.product(strategies, merging_args):
                data.append(measure(example, name, strategy, merge, args.repetitions))
                dframe = pd.DataFrame(data, index=mindex[:len(data)], columns=col_index)

                # Data is written after every measurement in case the job is killed
                dframe.to_csv(output_file)


    elif args.mode == "generate_code":

        linnea.config.set_verbosity(2)

        for strategy_str in strategy_strs:
            dir = os.path.join(linnea.config.results_path, args.experiment, "intensity", strategy_str)
            if not os.path.exists(dir):
                os.makedirs(dir)

        for example, name in job_examples:

            for strategy in strategies:
                generate(args.experiment, example, name, strategy)

            if args.reference:
                reference_code.generate_reference_code(name, example.eqns)
                operand_generation.generate_operand_generator(name, example.eqns)

                # runner should only include files that actually exists
                existing_algorithms = []
                for subdir_name, algorithm_name in [("constructive", "algorithm0c"), ("exhaustive", "algorithm0e")]:
                    file_path = os.path.join(linnea.config.output_code_path, name, Language.Julia.name, subdir_name, algorithm_name + ".jl")
                    if os.path.exists(file_path):
                        existing_algorithms.append((subdir_name, algorithm_name))
                
                runner.generate_runner(name, existing_algorithms)

    elif args.mode == "jobscripts":
        generate_scripts(args.experiment, len(examples))


if __name__ == "__main__":
    main()
