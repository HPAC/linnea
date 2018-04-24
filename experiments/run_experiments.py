import pandas as pd
import time
import numpy
import itertools
import sys
import math
import argparse

import linnea.config

linnea.config.set_output_path("~/linnea/output/")
linnea.config.set_verbosity(0)
linnea.config.init()

from linnea.config import Strategy
import linnea.examples.examples
import linnea.examples.lamp_paper.examples as lamp_paper
from linnea.derivation.graph.constructive import DerivationGraph as CDGraph
from linnea.derivation.graph.exhaustive import DerivationGraph as EDGraph

def measure(example, strategy, merge, job_index, reps=10):

    times = []

    if strategy is Strategy.constructive:
        DerivationGraph = CDGraph
    elif strategy is Strategy.exhaustive:
        DerivationGraph = EDGraph
    else:
        raise NotImplementedError()   

    for i in range(reps):
        linnea.config.clear_all()
        graph = DerivationGraph(example.eqns)
        t_start = time.perf_counter()
        trace = graph.derivation(
                            solution_nodes_limit=100,
                            iteration_limit=15,
                            merging=merge,
                            dead_ends=True)
        graph.write_output(code=True,
                           pseudocode=False,
                           output_name="tmp{}".format(job_index),
                           operand_generator=False,
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

def main():

    parser = argparse.ArgumentParser(prog="run_exp")
    parser.add_argument("job_index", help="Job index.", type=int)
    parser.add_argument("repetitions", help="Number of repetitions.", type=int, default=10)
    parser.add_argument("-c", "--constructive", help="Use constructive strategy.")
    parser.add_argument("-e", "--exhaustive", help="Use exhaustive strategy.")
    parser.add_argument("-n", "--no_merging", help="Also use no merging.")

    job_index = parser.job_index-1
    reps = parser.repetitions

    strategies = []
    if parser.constructive and parser.exhaustive:
        strategies = [Strategy.exhaustive, Strategy.constructive]
    elif parser.constructive:
        strategies = [Strategy.constructive]
    elif parser.exhaustive:
        strategies = [Strategy.exhaustive]
    else:
        return

    merging_args = [True]
    merging_labels = ["merging"]
    if parser.no_merging:
        merging_args.append(False)
        merging_labels.append("no_merging")


    # TODO doesn't work this way. Technically, init has to be called anew for each repetition.
    # Obviously, this is a problem if arguments are passed.
    # Reason: stuff like special properties. Properties in general?
    experiments = [
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

    data = []
    for strategy, merge in itertools.product(strategies, merging_args):
        data.append(measure(experiments[job_index], strategy, merge, job_index, reps))

    mindex = pd.MultiIndex.from_product([[type(experiments[job_index]).__name__], [strategy.name for strategy in strategies], merging_labels], names=["example", "strategy", "merging"])
    col_index = pd.Index(["mean", "std", "min", "max", "nodes"])

    dframe = pd.DataFrame(data, index=mindex, columns=col_index)
    dframe.to_pickle("linnea_generation{}.pkl".format(job_index))
    # print(dframe)


if __name__ == "__main__":
    main()
