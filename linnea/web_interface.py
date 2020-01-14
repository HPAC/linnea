
from .frontend.utils import parse_input
from .derivation.graph.derivation import DerivationGraph

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
