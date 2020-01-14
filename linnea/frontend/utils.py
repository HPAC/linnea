from tatsu.model import ModelBuilderSemantics

from .AST_translation import LinneaWalker
from .parser import LinneaParser

def parse_input(input):
    """Parses the input.

    Args:
        input (string): Description of the input in the Linnea language.

    Returns:
        Equations: The input equations.
    """
    parser = LinneaParser(semantics=ModelBuilderSemantics())
    ast = parser.parse(input, rule_name = "model")

    walker = LinneaWalker()
    walker.walk(ast)
    return walker.equations