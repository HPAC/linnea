
import argparse

parser = argparse.ArgumentParser(prog="linnea")
parser.add_argument("input", help="Relative path to the definition of the input.")
args = parser.parse_args()

with open(args.input, "r" ) as input_file:

    from . import config

    config.set_language(config.Language.Julia)
    config.set_data_type(config.JuliaDataType.Float64)
    config.init()

    from .derivation.graph.constructive import DerivationGraph
    # from .derivation.graph.exhaustive import DerivationGraph

    from . import examples

    from grako.model import ModelBuilderSemantics
    from .frontend.AST_translation import LinneaWalker
    from .frontend.parser import LinneaParser

    parser = LinneaParser(semantics=ModelBuilderSemantics())
    ast = parser.parse(input_file.read(), rule_name = "model")

    walker = LinneaWalker()
    walker.walk(ast)

    graph = DerivationGraph(walker.equations)
    trace = graph.derivation(10, 6)
    # print(":".join(str(t) for t in trace))

    graph.to_dot_file()
    graph.write_output(code=True,
                       pseudocode=True,
                       # output_name=args.input.split(".")[0],
                       output_name="tmp",
                       operand_generator=True,
                       max_algorithms=100,
                       graph=True)