
import argparse

parser = argparse.ArgumentParser(prog="linnea")
group = parser.add_mutually_exclusive_group()
parser.add_argument("input", help="input file")
parser.add_argument("output_path", nargs="?", help="relative path to the output directory; defaults to the current directory")
parser.add_argument("-d", "--data-type", choices=["Float32", "Float64"], help="data type used in the generated code")
parser.add_argument("-e", "--experiments", action="store_const", const=True, help="generate code for experiments")
parser.add_argument("-g", "--graph", action="store_const", const=True, help="write graph to file")
parser.add_argument("--max-algorithms", type=int, help="maximum number of generated algorithms")
parser.add_argument("--no-code", action="store_const", const=True, help="disable code generation")
parser.add_argument("--no-merging", action="store_const", const=True, help="disable merging branches in the graph")
parser.add_argument("-o", "--output", help="name of the output")
parser.add_argument("--pseudocode", action="store_const", const=True, help="generate pseudocode")
group.add_argument("--silent", action="store_const", const=True, help="suppress non-error messages")
parser.add_argument("-s", "--strategy", choices=["constructive", "exhaustive"], help="specify derivation strategy")
group.add_argument("-v", "--verbose", action="count", help="increase verbosity")
parser.add_argument("--example", help=argparse.SUPPRESS)
args = parser.parse_args()


from . import config

# config is used to make sure that the values from config.json are also considered.

if args.output_path is not None:
    config.set_output_path(args.output_path)
if args.data_type is not None:
    config.set_data_type(config.JuliaDataType[args.data_type])
if args.experiments is not None:
    config.set_generate_experiments(args.experiments)
if args.graph is not None:
    config.set_generate_graph(args.graph)
if args.max_algorithms is not None:
    config.set_max_algorithms(args.max_algorithms)
if args.code is not None:
    config.set_generate_code(args.code)
if args.merging is not None:
    config.set_merging_branches(args.merging)
if args.output is not None:
    config.set_output_name(args.output)
else:
    config.set_output_name(os.path.splitext(os.path.basename(args.input))[0])
if args.output_path is not None:
    config.set_output_path(args.output_path)
if args.pseudocode is not None:
    config.set_generate_pseudocode(args.pseudocode)
if args.silent is not None:
    config.set_verbosity(0)
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
        from grako.model import ModelBuilderSemantics
        from .frontend.AST_translation import LinneaWalker
        from .frontend.parser import LinneaParser

        parser = LinneaParser(semantics=ModelBuilderSemantics())
        ast = parser.parse(input_file.read(), rule_name = "model")

        walker = LinneaWalker()
        walker.walk(ast)
        equations = walker.equations
else:
    from . import examples
    # TODO retrieve example

graph = DerivationGraph(equations)
trace = graph.derivation(10, 6, merging=config.merging_branches)
if config.verbosity >= 2:
    print(":".join(str(t) for t in trace))

graph.write_output(code=config.generate_code,
                   pseudocode=config.generate_pseudocode,
                   output_name=config.output_name,
                   operand_generator=config.generate_experiments,
                   max_algorithms=config.max_algorithms,
                   graph=config.generate_graph)