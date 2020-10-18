from .frontend.utils import parse_input
from .frontend.export import export_expression
from .derivation.graph.derivation import DerivationGraph
from . import utils

from linnea.algebra.validity import ExpressionException, InvalidExpression, \
                                    SizeMismatch, ConflictingProperties

import linnea.config
import linnea.algebra.expression as ae

import json
import re

class SyntaxError(Exception):
    pass

class GenerationError(Exception):
    pass

variable_regex = re.compile("([a-zA-Z_][a-zA-Z0-9_]*)\s*=\s*([0-9]+)")
matrix_reges = re.compile("(IdentityMatrix|ZeroMatrix)\s*([a-zA-Z_][a-zA-Z0-9_]*)\s*\(\s*([a-zA-Z_][a-zA-Z0-9_]*)\s*,\s*([a-zA-Z_][a-zA-Z0-9_]*)\s*\)")

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

    try:
        equations = parse_input(input)
    except:
        raise SyntaxError()
    
    try:
        graph = DerivationGraph(equations)
        graph.derivation(time_limit=time_limit, merging=True, pruning_factor=1.)
        output = graph.optimal_algorithm_to_str()
    except ConflictingProperties as e:
        raise e
    except ExpressionException as e:
        raise e.replace_expressions()
    except:
        raise GenerationError("An error occurred during the algorithm generation.")

    if not output:
        for equation in equations:
            for expr, _ in equation.rhs.preorder_iter():
                if isinstance(expr, ae.Vector) and expr.rows == 1:
                    raise GenerationError("Row vectors are not fully supported yet.")
        raise GenerationError("An error occurred during the algorithm generation.")

    return output


def dependent_dimensions(input):
    """Computes dependent dimensions.

    The dependent dimensions are all sets of dimensions that have to be the same
    for the input equations to be valid. The output is a JSON string, and
    dependent dimensions are represented as nested arrays. The innermost arrays
    represent the dimensions; they contain two elements: The first element is
    the name of the operand as a string, the second element is an integer; 0
    stands for rows, 1 for columns. All dimensions that have to be the same are
    in the same array.

    For the input, the custom input language of Linnea has to be used.
    
    Args:
        input (str): Description of the input.

    Returns:
        string: A JSON string of nested arrays.
    """
    try:
        equations = parse_input(input)
    except:
        raise SyntaxError()

    try:
        equations.check_validity(dependent_dimensions=True)
    except ExpressionException as e:
        raise e.replace_expressions()

    # The internal names of constant matrices are replaced with their names in
    # the input.
    variables = dict()    
    for name, val in variable_regex.findall(input):
        variables[name] = val

    renaming = dict()
    for type, name, rows, cols in matrix_reges.findall(input):
        if type == "IdentityMatrix":
            prefix = "I"
        elif type == "ZeroMatrix":
            prefix = "0"
        internal_name = "{}({}, {})".format(prefix, variables[rows], variables[cols])
        renaming[internal_name] = name

    dimensions = []
    for l in utils.dependent_dimensions(equations):
        new_l = []
        for name, dim in l:
            name = renaming.get(name, name)
            new_l.append([name, dim])
        dimensions.append(new_l)

    return json.dumps(dimensions)


