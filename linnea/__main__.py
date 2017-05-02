
import argparse

parser = argparse.ArgumentParser(prog="clak2")
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

    from .frontend.AST_translation import ASTTranslator
    from .frontend.parser import LinneaParser

    parser = LinneaParser()
    ast = parser.parse(input_file.read(), rule_name = "model")

    ast_translator = ASTTranslator(ast)

    graph = DerivationGraph(ast_translator.equations)
    trace = graph.derivation(10, 6)
    # print(":".join(str(t) for t in trace))

    graph.to_dot_file()
    graph.algorithms_to_files(pseudocode=True, experiment=args.input.split(".")[0])