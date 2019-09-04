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

from linnea.config import Language
import linnea.examples.examples
import linnea.examples.application as application
from linnea.derivation.graph.derivation import DerivationGraph
from linnea.code_generation.experiments.utils import generate_experiment_code

from random_expressions import generate_equation
from jobscripts import generate_scripts

def measure(experiment, example, name, merging):

    linnea.config.clear_all()
    if hasattr(example, "init"):
        # calls initialization that have to be done before each repetition
        example.init()

    graph = DerivationGraph(example.eqns)
    trace = graph.derivation(
                        time_limit=30*60,
                        merging=merging,
                        dead_ends=True)

    df_trace = pd.DataFrame(trace, columns=["time", "cost"])

    t_start = time.perf_counter()
    graph.write_output(code=True,
                       derivation=False,
                       output_name=name,
                       experiment_code=False,
                       algorithms_limit=1,
                       graph=False,
                       subdir_name="time_generation")
    t_end = time.perf_counter()

    mindex = pd.MultiIndex.from_tuples([(name, merging)], names=["example", "merging"])
    df_code_gen_time = pd.DataFrame([t_end-t_start], index=mindex, columns=["time"])
    
    if merging:
        subdir = "merging"
    else:
        subdir = "no_merging"

    file_path = os.path.join(linnea.config.results_path, experiment, "generation", subdir, name + "_trace.csv")
    df_trace.to_csv(file_path)
    if linnea.config.verbosity >= 2:
        print("Generate trace file {}".format(file_path))

    file_path = os.path.join(linnea.config.results_path, experiment, "generation", subdir, name + "_code_gen_time.csv")
    df_code_gen_time.to_csv(file_path)
    if linnea.config.verbosity >= 2:
        print("Generate code gen time file {}".format(file_path))

    return

def generate(experiment, example, name, k_best=False):

    algorithm_name = "algorithm{}"

    linnea.config.clear_all()
    if hasattr(example, "init"):
        # calls initialization that have to be done before each repetition
        example.init()

    print(example.eqns)

    if k_best:
        # time_limit = ???
        algorithms_limit = 100
        pruning_factor = 1.5
        intensity_dir = "intensity_k_best"
        subdir_name_experiments = "k_best"
    else:
        algorithms_limit = 1
        pruning_factor = 1.0
        intensity_dir = "intensity"
        subdir_name_experiments = "experiments"

    graph = DerivationGraph(example.eqns)
    trace = graph.derivation(
                        time_limit=30*60,
                        merging=True,
                        dead_ends=True,
                        pruning_factor=pruning_factor)

    k = graph.write_output(
                        code=False,
                        derivation=True,
                        output_name=name,
                        experiment_code=True,
                        algorithms_limit=algorithms_limit,
                        graph=False,
                        subdir_name_experiments=subdir_name_experiments,
                        algorithm_name=algorithm_name,
                        no_duplicates=True)

    vals = []
    data = example.eqns.get_data()
    for _, cost in graph.k_shortest_paths(k):
        vals.append([data, cost, cost/data])

    mindex = pd.MultiIndex.from_product([[name], [algorithm_name.format(i) for i in range(k)]], names=("example", "implementation"))
    df = pd.DataFrame(vals, index=mindex, columns=["data", "cost", "intensity"])

    file_path = os.path.join(linnea.config.results_path, experiment, intensity_dir, name + "_intensity.csv")
    df.to_csv(file_path)
    if linnea.config.verbosity >= 2:
        print("Generate intensity file {}".format(file_path))

def main():

    parser = argparse.ArgumentParser(prog="experiments")
    parser.add_argument("mode", choices=["time_generation", "generate_code", "jobscripts"])
    parser.add_argument("experiment", choices=["random", "application"])
    parser.add_argument("-m", "--merging", choices=["true", "false", "both"], default="true")
    parser.add_argument("-j", "--jobindex", help="Job index.", type=int, default=0)
    parser.add_argument("-f", "--reference", action="store_true", help="Only generate reference code.")
    parser.add_argument("--k-best", action="store_true", help="Generate code for k best experiment.")
    parser.add_argument("-l", "--config", type=str, default=None, help="Specify configuration file.")
    args = parser.parse_args()

    application_examples = [
        application.Example01(),
        application.Example02(),
        application.Example03(),
        application.Example04(),
        application.Example05(),
        application.Example06(),
        application.Example07(),
        application.Example08(),
        application.Example09(),
        application.Example10(),
        application.Example11(),
        application.Example12(),
        application.Example13(),
        application.Example14(),
        application.Example15(),
        application.Example16(),
        application.Example17(),
        application.Example18(),
        application.Example19(),
        application.Example20(),
        application.Example21(),
        application.Example22(),
        application.Example23(),
        application.Example24(),
        application.Example25(),
    ]

    ExampleContainer = collections.namedtuple("ExampleContainer", ["eqns"])

    random.seed(0)
    rand_exprs = [ExampleContainer(generate_equation(random.randint(4, 7))) for _ in range(100)]

    if args.config:
        linnea.config.load_config(config_file=args.config)

    if args.experiment == "random":
        examples = rand_exprs
    elif args.experiment == "application":
        examples = application_examples
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


    num_threads = [1, 24]

    if args.mode == "time_generation":

        linnea.config.set_verbosity(0)

        # TODO this is not necessary because it also happens in the jobscript. It might be possible to simplify the jobscript too.
        # All paths are known here, there might be no need to go anywhere in the jobscript.
        for subdir in ["merging", "no_merging"]:
            dir = os.path.join(linnea.config.results_path, args.experiment, "generation", subdir)
            if not os.path.exists(dir):
                os.makedirs(dir, exist_ok=True) # exist_ok=True avoids errors when running experiments in parallel

        merging_args = []
        if args.merging == "true" or args.merging == "both":
            merging_args.append(True)
        elif args.merging == "false" or args.merging == "both":
            merging_args.append(False)

        data = []
        for example, name in job_examples:
            for merge in merging_args:
                measure(args.experiment, example, name, merge)


    elif args.mode == "generate_code":

        linnea.config.set_verbosity(2)

        for subdir in ["trace", "intensity", "intensity_k_best"]:
            dir = os.path.join(linnea.config.results_path, args.experiment, subdir)
            if not os.path.exists(dir):
                os.makedirs(dir, exist_ok=True) # exist_ok=True avoids errors when running experiments in parallel

        for example, name in job_examples:

            if not args.reference:
                generate(args.experiment, example, name, args.k_best)
            else:
                algorithm_name = "algorithm{}"

                # determining k for k best experiment
                # using the upper limit for k is sufficient because we test if files exist
                k = 0
                while True:
                    file_path = os.path.join(linnea.config.output_code_path, name, Language.Julia.name, "k_best", algorithm_name.format(k) + ".jl")
                    if os.path.exists(file_path):
                        k += 1
                    else:
                        break

                generate_experiment_code(name, example.eqns, algorithm_name, k, num_threads)

    elif args.mode == "jobscripts":
        generate_scripts(args.experiment, len(examples), num_threads)


if __name__ == "__main__":
    main()
