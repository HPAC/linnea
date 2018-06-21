import pandas as pd
import time
import numpy
import itertools
import sys
import math
import argparse

import linnea.config

linnea.config.set_output_path("~/linnea/output/")
linnea.config.set_verbosity(2)
linnea.config.init()

from linnea.config import Strategy
import linnea.examples.examples
import linnea.examples.lamp_paper.examples as lamp_paper
from linnea.derivation.graph.constructive import DerivationGraph as CDGraph
from linnea.derivation.graph.exhaustive import DerivationGraph as EDGraph
from linnea.code_generation.experiments import operand_generation, runner, reference_code
from linnea.code_generation import utils as cgu

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
                           pseudocode=False,
                           output_name=name,
                           experiment_code=False,
                           algorithms_limit=1,
                           graph=False)
        t_end = time.perf_counter()
        times.append(t_end-t_start)
    if graph.terminal_nodes():
        data = [numpy.mean(times), numpy.std(times), numpy.min(times), numpy.max(times)]
        data.append(graph.nodes_count())
    else:
        data = [math.nan]*5
    return data

def generate(example, name, strategy):

    if strategy is Strategy.constructive:
        DerivationGraph = CDGraph
        algorithm_name = "algorithm{}c"
    elif strategy is Strategy.exhaustive:
        DerivationGraph = EDGraph
        algorithm_name = "algorithm{}e"
    else:
        raise NotImplementedError()

    graph = DerivationGraph(example.eqns)
    trace = graph.derivation(
                        solution_nodes_limit=100,
                        iteration_limit=15,
                        merging=True,
                        dead_ends=True)
    graph.write_output(code=True,
                       pseudocode=False,
                       output_name=name,
                       experiment_code=False,
                       algorithms_limit=1,
                       graph=False,
                       subdir_name=strategy.name,
                       algorithm_name=algorithm_name)

def generate_name(index):
    return "lamp_example{}".format(index)

def main():

    parser = argparse.ArgumentParser(prog="experiments")
    parser.add_argument("mode", choices=["time_generation", "generate_code"])
    parser.add_argument("-j", "--jobindex", dest="jobindex", help="Job index.", type=int, default=0)
    parser.add_argument("-r", "--repetitions", help="Number of repetitions.", type=int)
    parser.add_argument("-c", "--constructive", action="store_true", help="Use constructive strategy.")
    parser.add_argument("-e", "--exhaustive", action="store_true", help="Use exhaustive strategy.")
    parser.add_argument("-n", "--no-merging", dest="no_merging", action="store_true", help="Also use no merging.")
    args = parser.parse_args()

    lamp_examples = [
        lamp_paper.LeastSquares_7_1_1(),
        lamp_paper.LMMSE_7_1_2(),
        lamp_paper.Generalized_LeastSquares_7_1_3(),
        lamp_paper.Optimization_Problem_7_1_4(),
        lamp_paper.Signal_Processing_7_1_5(),
        lamp_paper.Lower_Triangular_Inversion_7_1_6(),
        lamp_paper.Local_Assimilation_Kalman_7_1_7(),
        lamp_paper.EnsembleKalmanFilter_7_1_9_1(),
        lamp_paper.EnsembleKalmanFilter_7_1_9_2(),
        lamp_paper.SPA_7_1_12(),
        lamp_paper.SPA_7_1_12(q=2),
        lamp_paper.ImageRestoration_7_1_13_1(),
        lamp_paper.ImageRestoration_7_1_13_2(single=False),
        lamp_paper.ImageRestoration_7_1_13_2(single=True),
        lamp_paper.Tikhonov_7_1_14(),
        lamp_paper.CDMA_7_1_15(630, 300, 50),
        lamp_paper.Common_Subexpr_7_2_1(),
        lamp_paper.Common_Subexpr_7_2_2(),
        lamp_paper.Common_Subexpr_7_2_3(),
        lamp_paper.Overlap_Common_Subexpr_7_2_4(),
        lamp_paper.Rewrite_Distributivity_7_2_5_1(),
        lamp_paper.Rewrite_Distributivity_7_2_5_2(),
        lamp_paper.Rewrite_Distributivity_7_2_5_3(),
        lamp_paper.Rewrite_Distributivity_7_2_5_4(),
        lamp_paper.Matrix_Chain_7_2_6(),
        lamp_paper.Matrix_Chain_7_2_7(),
        lamp_paper.Properties_7_2_8(),
        lamp_paper.Transposed_Kernel_7_2_9(),
        lamp_paper.Transposed_Kernel_7_2_10(),
        lamp_paper.Simplification_7_2_11(),
        lamp_paper.Simplification_7_2_12(),
    ]

    # TODO when using different sets of examples, use that for output name
    # also add init for new examples

    # TODO write parameters (repetitions) to file
    # TODO write description of experiment (name, operand sizes, equations) to file

    if args.jobindex == 0:
        examples = enumerate(lamp_examples, 1)
    else:
        examples = [(args.jobindex, lamp_examples[args.jobindex-1])]

    strategies = []
    if args.constructive:
        strategies.append(Strategy.constructive)
    if args.exhaustive:
        strategies.append(Strategy.exhaustive)
    if not strategies:
        return

    if args.mode == "time_generation":

        merging_args = [True]
        merging_labels = ["merging"]
        if args.no_merging:
            merging_args.append(False)
            merging_labels.append("no_merging")

        data = []
        for idx, example in examples:
            for strategy, merge in itertools.product(strategies, merging_args):
                name = generate_name(idx, strategy)
                data.append(measure(example, name, strategy, merge, args.repetitions))

        mindex = pd.MultiIndex.from_product([[type(exp).__name__ for exp in examples], [strategy.name for strategy in strategies], merging_labels], names=["example", "strategy", "merging"])
        col_index = pd.Index(["mean", "std", "min", "max", "nodes"])

        dframe = pd.DataFrame(data, index=mindex, columns=col_index)

        if args.jobindex == 0:
            dframe.to_csv("linnea_generation.csv")
        else:  
            dframe.to_csv("linnea_generation{}.csv".format(args.jobindex))
        # print(dframe)

    elif args.mode == "generate_code":
        for idx, example in examples:
            name = generate_name(idx)
            for strategy in strategies:
                generate(example, name, strategy)

            reference_code.generate_reference_code(name, example.eqns)
            operand_generation.generate_operand_generator(name, example.eqns)

            cgu.remove_files(os.path.join(config.output_path, name, config.language.name, "constructive"))
            cgu.remove_files(os.path.join(config.output_path, name, config.language.name, "exhaustive"))

            algorithms = [("constructive", "algorithm0c"), ("exhaustive", "algorithm0e")]
            runner.generate_runner(name, algorithms)




if __name__ == "__main__":
    main()
