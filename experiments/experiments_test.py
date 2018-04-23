import pandas as pd
import time
import numpy
import itertools
import sys

import linnea.config

linnea.config.set_output_path("~/linnea/output/")
linnea.config.set_verbosity(0)
linnea.config.init()

from linnea.config import Strategy
import linnea.examples.examples
import linnea.examples.lamp_paper.examples as lamp_paper
from linnea.derivation.graph.constructive import DerivationGraph as CDGraph
from linnea.derivation.graph.exhaustive import DerivationGraph as EDGraph

def measure(example, strategy, merge, reps=10):

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
                           output_name="tmp",
                           operand_generator=False,
                           algorithms_limit=1,
                           graph=False)
        t_end = time.perf_counter()
        times.append(t_end-t_start)
    data = [numpy.mean(times), numpy.std(times), numpy.min(times), numpy.max(times)]
    data.append(graph.nodes_count())
    return data

if __name__ == "__main__":

    if len(sys.argv) > 1:
        reps = int(sys.argv[1])

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
    for experiment in experiments:
        for strategy, merge in itertools.product([Strategy.exhaustive, Strategy.constructive], [True, False]):
            data.append(measure(experiment, strategy, merge, reps))

    mindex = pd.MultiIndex.from_product([[type(ex).__name__ for ex in experiments], ["exhaustive", "constructive"], ["merging", "no_merging"]], names=["example", "strategy", "merging"])
    col_index = pd.Index(["mean", "std", "min", "max", "nodes"])

    dframe = pd.DataFrame(data, index=mindex, columns=col_index)
    dframe.to_csv("linnea_generation.csv")
    # print(dframe)
