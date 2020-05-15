from .frontend.utils import parse_input
from .derivation.graph.derivation import DerivationGraph
from . import utils

import linnea.config

def run_linnea(input, time_limit=10):
    """Run Linnea code generation.

    This function runs the code generation of Linnea on the input and returns
    the code that implements the optimal algorithm as a string.

    For the input, the custom input language of Linnea has to be used.
    
    Args:
        input (str): Description of the input.
        time_limit (int, optional): Time limit for the generation in seconds.

    Returns:
        str: The generated code.
    """
    
    linnea.config.set_verbosity(0)

    equations = parse_input(input)

    graph = DerivationGraph(equations)
    graph.derivation(time_limit=time_limit, merging=True, pruning_factor=1.)

    return graph.optimal_algorithm_to_str()


def dependent_dimensions(input):
    """Computes dependent dimensions.

    The dependent dimensions are all sets of dimensions that have to be the same
    for the input equations to be valid. Dimensions are represented as tuples of
    two elements: The first element is the name of the operand, the second
    element is an integer; 0 stands for rows, 1 for columns.

    For the input, the custom input language of Linnea has to be used.
    
    Args:
        input (str): Description of the input.

    Returns:
        list: A list of sets. All dimensions in one set have to be the same.
    """
    equations = parse_input(input)
    return utils.dependent_dimensions(equations)