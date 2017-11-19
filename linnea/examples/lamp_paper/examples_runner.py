
from timeit import default_timer as timer

import numpy, csv, math

from linnea import config
from linnea.code_generation import experiments as cge
from linnea.derivation.graph.constructive import DerivationGraph

from linnea.examples.examples import Example126
from linnea.examples.examples import Example127
from linnea.examples.lamp_paper import examples as paper

def generate(eq, append, merging_, algos = 100, iters = 10,reps=10):
    times = []
    #graph = DerivationGraph(eq)
    # it used to be (10,6)
    #trace = graph.derivation(algos, iters, verbose=False, merging=merging_)
    for j in range(reps):
        start = timer()
        graph = DerivationGraph(eq)
        trace = graph.derivation(algos, iters, verbose=False, merging=merging_)
        end = timer()
        times.append(end - start)
    time_measurements = [numpy.mean(times), numpy.std(times), numpy.min(times), numpy.max(times)]
    graph = DerivationGraph(eq)
    trace = graph.derivation(algos, iters, verbose=False, merging=merging_)
    graph.write_output(code=True,
                       #pseudocode=True,
                       # output_name=args.input.split(".")[0],
                       output_name="example_{}_{}_{}".format(type(equations).__name__, append, "merging" if merging_  else "nomerging"),
                       operand_generator=True,
                       max_algorithms=100)
                       #graph=True)
    time_measurements.append(graph.nodes_count())
    return time_measurements
    

if __name__ == "__main__":

    import sys
    algorithm = int(sys.argv[1]) - 1
    config.set_language(config.Language.Julia)
    config.set_data_type(config.JuliaDataType.Float64)
    config.init()


    dirs = []
    eqs = []
    #eqs.append(paper.LeastSquares_7_1_1())
    eqs.append(paper.LMMSE_7_1_2())
    #eqs.append(paper.Generalized_LeastSquares_7_1_3())
    #eqs.append(paper.Optimization_Problem_7_1_4())
    ## very very long
    ##eqs.append(paper.Signal_Processing_7_1_5())
    #eqs.append(paper.Lower_Triangular_Inversion_7_1_6())
    # bug
    #eqs.append(paper.Local_Assimilation_Kalmar_7_1_7())
    eqs.append(paper.EnsembleKalmarFilter_7_1_9_1())
    eqs.append(paper.EnsembleKalmarFilter_7_1_9_2())
    #eqs.append(paper.SPA_7_1_12())
    #eqs.append(paper.SPA_7_1_12(q=2))
    #  IndexError: pop index out of range
    #eqs.append(paper.ImageRestoration_7_1_13_1())
    #eqs.append(paper.ImageRestoration_7_1_13_2(single=False))
    # bug 
    #eqs.append(paper.ImageRestoration_7_1_13_2(single=True))
    eqs.append(paper.Tikhonov_7_1_14())
    eqs.append(paper.CDMA_7_1_15(630, 300, 50))
    eqs.append(paper.Common_Subexpr_7_2_1())
    eqs.append(paper.Common_Subexpr_7_2_2())
    eqs.append(paper.Common_Subexpr_7_2_3())
    eqs.append(paper.Overlap_Common_Subexpr_7_2_4())
    eqs.append(paper.Rewrite_Distributivity_7_2_5_1())
    eqs.append(paper.Rewrite_Distributivity_7_2_5_2())
    eqs.append(paper.Rewrite_Distributivity_7_2_5_3())
    eqs.append(paper.Rewrite_Distributivity_7_2_5_4())
    eqs.append(paper.Matrix_Chain_7_2_6())
    eqs.append(paper.Matrix_Chain_7_2_7())
    eqs.append(paper.Properties_7_2_8())
    eqs.append(paper.Transposed_Kernel_7_2_9())
    eqs.append(paper.Transposed_Kernel_7_2_10())
    eqs.append(paper.Simplification_7_2_11())
    eqs.append(paper.Simplification_7_2_12())
    equations = eqs[algorithm] 
        #print(i)
        #equations = next(gen) 
        #equations = eqs[i]
    print(equations.eqns)
    measurements = []
    titles = []
    equation_name = type(equations).__name__
    
    from linnea.derivation.graph.constructive import DerivationGraph
    dirs.append("example_{}_{}_{}".format(type(equations).__name__, "constructive","merging"))
    dirs.append("example_{}_{}_{}".format(type(equations).__name__, "constructive","nomerging"))
    measurements.append( generate(equations.eqns, "constructive", False) )
    titles.append("constructive_nonmerging")
    measurements.append( generate(equations.eqns,"constructive", True) )
    titles.append("constructive_merging")

    from linnea.derivation.graph.exhaustive import DerivationGraph
    dirs.append("example_{}_{}_{}".format(type(equations).__name__, "exhaustive", "merging"))
    dirs.append("example_{}_{}_{}".format(type(equations).__name__, "exhaustive", "nomerging"))
    measurements.append( generate(equations.eqns, "exhaustive", False, reps=5) )
    titles.append("exhaustive_nonmerging")
    measurements.append( generate(equations.eqns, "exhaustive", True, reps=5) )
    titles.append("exhaustive_merging")
    #FIXME: independent from language
    with open('example_{}_time_measurements'.format(equation_name), 'w') as f:
        writer = csv.writer(f,delimiter="\t")
        writer.writerow(["algorithm","Time","StdDev","Min","Max","Nodes_count"])
        for row_title, data_row in zip(numpy.array(titles), numpy.array(measurements)):
             writer.writerow([row_title] + data_row.tolist())
    #cge.generate_cmake_script(dirs)
