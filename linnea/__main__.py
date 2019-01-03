
import argparse
import os
import operator
import logging

logging.basicConfig(level=logging.INFO,
                    format='%(name)-2s: %(levelname)-2s %(message)s')

from . import config

class ExampleDoesNotExist(Exception):
    pass

def main():

    parser = argparse.ArgumentParser(prog="linnea")
    group = parser.add_mutually_exclusive_group()
    parser.add_argument("input", help="input file")
    parser.add_argument("output_code_path", nargs="?", help="relative path to the output directory; defaults to the current directory")
    parser.add_argument("--algorithms-limit", type=int, help="maximum number of generated algorithms")
    parser.add_argument("-d", "--data-type", choices=["Float32", "Float64"], help="data type used in the generated code")
    parser.add_argument("--no-dead-ends", dest="dead_ends", action="store_const", const=True, help="disable dead end detection")
    parser.add_argument("-e", "--experiments", action="store_const", const=True, help="generate code for experiments")
    parser.add_argument("-g", "--graph", action="store_const", const=True, help="write graph to file")
    parser.add_argument("--graph-style", choices=["full", "simple", "minimal"], help="specify graph style")
    parser.add_argument("--iteration-limit", help="iteration limit")
    parser.add_argument("--no-code", dest="code", action="store_const", const=True, help="disable code generation")
    parser.add_argument("--no-merging", dest="merging", action="store_const", const=True, help="disable merging branches in the graph")
    parser.add_argument("-o", "--output", help="name of the output")
    parser.add_argument("--derivation", action="store_const", const=True, help="generate description of the derivation")
    group.add_argument("--silent", action="store_const", const=True, help="suppress non-error messages")
    parser.add_argument("--solution-nodes-limit", help="limit for the number of solution nodes")
    parser.add_argument("-s", "--strategy", choices=["constructive", "exhaustive"], help="specify derivation strategy")
    group.add_argument("-v", "--verbose", action="count", help="increase verbosity")
    parser.add_argument("--example", help=argparse.SUPPRESS)
    args = parser.parse_args()

    # config is used to make sure that the values from config.json are also considered.
    # print(parser.args)

    if args.output_code_path is not None:
        config.set_output_code_path(args.output_code_path)
    if args.algorithms_limit is not None:
        config.set_algorithms_limit(args.algorithms_limit)
    if args.data_type is not None:
        config.set_data_type(config.JuliaDataType[args.data_type])
    if args.experiments is not None:
        config.set_generate_experiments(args.experiments)
    if args.dead_ends is not None:
        config.set_dead_ends(args.dead_ends)
    if args.graph is not None:
        config.set_generate_graph(args.graph)
    if args.graph_style is not None:
        config.set_graph_style(config.GraphStyle[args.graph_style])
    if args.iteration_limit is not None:
        config.set_iteration_limit(args.iteration_limit)
    if args.code is not None:
        config.set_generate_code(args.code)
    if args.merging is not None:
        config.set_merging_branches(args.merging)
    if args.output is not None:
        config.set_output_name(args.output)
    else:
        config.set_output_name(os.path.splitext(os.path.basename(args.input))[0])
    if args.derivation is not None:
        config.set_generate_derivation(args.derivation)
    if args.silent is not None:
        config.set_verbosity(0)
    if args.solution_nodes_limit is not None:
        config.set_solution_nodes_limit(args.solution_nodes_limit)
    if args.strategy is not None:
        config.set_strategy(config.Strategy[args.strategy])
    if args.verbose is not None:
        config.set_verbosity(args.verbose+1)

    config.init()

    if config.strategy is config.Strategy.constructive:
        from .derivation.graph.constructive import DerivationGraph  
    elif config.strategy is config.Strategy.exhaustive:
        from .derivation.graph.exhaustive import DerivationGraph


    if args.example is None:
        with open(args.input, "r" ) as input_file:
            from tatsu.model import ModelBuilderSemantics
            from .frontend.AST_translation import LinneaWalker
            from .frontend.parser import LinneaParser

            parser = LinneaParser(semantics=ModelBuilderSemantics())
            ast = parser.parse(input_file.read(), rule_name = "model")

            walker = LinneaWalker()
            walker.walk(ast)
            equations = walker.equations

            # print(equations)
            # quit()

    else:
        from .examples import examples

        example_name = args.example
        try:
            example_idx = int(args.example)
        except ValueError:
            pass
        else:
            example_name = "Example{:03}".format(example_idx)

        caller = operator.methodcaller(example_name)
        try:
            example = caller(examples)
        except AttributeError:
            raise ExampleDoesNotExist(example_name)

        config.set_output_name(example_name)
        equations = example.eqns

    graph = DerivationGraph(equations)
    trace = graph.derivation(
                        solution_nodes_limit=config.solution_nodes_limit,
                        iteration_limit=config.iteration_limit,
                        merging=config.merging_branches)
    if config.verbosity >= 2:
        print(":".join(str(t) for t in trace))

    graph.write_output(code=config.generate_code,
                       derivation=config.generate_derivation,
                       output_name=config.output_name,
                       experiment_code=config.generate_experiments,
                       algorithms_limit=config.algorithms_limit,
                       graph=config.generate_graph,
                       graph_style=config.graph_style)

if __name__ == '__main__':
    main()
