from ...algebra.expression import Symbol, Matrix, Scalar, Inverse, \
                                  Times, Plus, Operator, \
                                  Transpose, ConjugateTranspose, \
                                  InverseTranspose, InverseConjugate, \
                                  InverseConjugateTranspose

from ...algebra.transformations import transpose, invert
from ...algebra.representations import to_POS

from ...algebra import equations as aeq

from ...algebra.properties import Property as properties

from collections import namedtuple, deque

from enum import Enum, unique

import copy
import itertools
import operator

@unique
class ExpressionType(Enum):
    # An expression of the form Inverse(Symbol), where Symbol has a property
    # that makes factorizations unnecessary.
    simple_inverse_no_factor = 1

    # An expression of the form Inverse(Symbol), where Symbol has a property
    # that makes factorizations necessary.
    simple_inverse = 2

    # An expression of the form Inverse(Times(…)) or Inverse(Plus(…)), where
    # the operands of Times/Plus have properties that make factorizations
    # unnecessary.
    compound_inverse_no_factor = 3

    # An expression of the form Inverse(Times(…)) or Inverse(Plus(…)), where
    # the operands of have properties that make factorizations necessary.
    compound_inverse = 4

    # One of the following:
    # Plus(…) where all operands are either Symbol, Transpose(Symbol) or
    # Minus(Symbol).
    # Times(…) where all operands are either Symbol, Transpose(Symbol), 
    # Minus(Symbol), … or simple inverses that don't require factorizations.
    simple = 5

    none = 6

@unique
class OperationType(Enum):
    plus = 1
    times = 2
    none = 3

# @profile
def generate_variants(equations, eqn_idx=None):
    """Generates "product of sums" variants of equations.

    This function generates some "product of sums" variants of the given
    equations and returns all unique ones, including the one that was passed to
    this function as an argument.

    The optional argument eqn_idx can be used to specify that variants are
    generated for only one single equation.
    """
    # If there is no Plus in any equation, there will be no POS variant.
    plus_found = False
    for equation in equations:
        for node, pos in equation.preorder_iter():
            if isinstance(node, Plus):
                plus_found = True
                break
        if plus_found:
            break
    else:
        return [equations]


    if eqn_idx is None:
        POSL_equations = aeq.Equations(*[to_POS(equation, "l") for equation in equations])
        POSR_equations = aeq.Equations(*[to_POS(equation, "r") for equation in equations])
    else:
        POSL_equations = copy.deepcopy(equations)
        POSR_equations = copy.deepcopy(equations)
        POSL_equations[eqn_idx] = to_POS(POSL_equations[eqn_idx], "l")
        POSR_equations[eqn_idx] = to_POS(POSR_equations[eqn_idx], "r")  

    use_POSL = True
    use_POSR = True

    if POSL_equations == equations:
        use_POSL = False
    if POSR_equations == equations:
        use_POSR = False
    if use_POSR and use_POSL and POSR_equations == POSL_equations:
        # If POSR and POSL equations are the same, only one is used.
        use_POSR = False

    variants = [equations]
    if use_POSR:
        variants.append(POSR_equations)
    if use_POSL:
        variants.append(POSL_equations)

    # print(variants)
    # return [equations]
    return variants


def process_next(equation):
    """Finds a subexpression to process next in the derivation.

    This function takes an equation as input. It is assumed that the equation
    has the form <symbol> = <expression>.

    Returns a tuple with
    - the path to the expression to process next
    - the type of that expression (as ExpressionType)
    - the operation type of that expression, if applicable (as OperationType)

    Subexpressions are selected as follows:
    1. The first inverse (according to postorder) that requires further
       processing.
    2. If there is no such inverese, the first (according to postorder) simple
       plus of simple times.
    """

    # It is necessary to search for inverses first because it's possible that
    # simple plus or times is inside an inverse, so factorizations have to be
    # used.
    for pos in inverse_positions(equation.rhs, [1]):
        node = equation[pos]
        expr_type, op_type = inverse_type(node)
        if expr_type == ExpressionType.simple_inverse_no_factor:
            continue
        else:
            return (pos, expr_type, op_type)


    # for node, pos in equation.rhs.iterate_postorder_with_positions([1]):
    for node, pos in process_next_generator(equation.rhs, [1]):
        # node = equation[pos]
        # print("here", node, node.size)
        # print(list((n, n.size) for n in node.iterate_preorder()))
        # if is_scalar(node):
            # print("here")
        if isinstance(node, Plus) and is_simple_plus(node):
            return (pos, ExpressionType.simple, OperationType.plus)
        elif isinstance(node, Times) and is_simple_times(node):
            return (pos, ExpressionType.simple, OperationType.times)

    # [1] = position of right-hand side
    return ([1], ExpressionType.none, OperationType.none)

def process_next_simple(expression, position=[]):
    """Finds a subexpression to process next in the derivation.

    This is a simplified verision of process_next(). The difference is that this
    function does not take inverses into account. Thus, is can be used for
    expression where it is already known that factorizations have to be applied
    (like compund inverses), but it still has to be determined what to do in
    addition to that (solving sums or products, or applying unary operators).

    As input, this function takes an expression. If this expression is a
    subexpression of a another expression, the optional argument position can be
    used to ensure that process_next_simple returns the correct path to the
    subexpression to process next.

    Returns a tuple with
    - the path to the expression to process next
    - the operation type of that expression, if applicable (as OperationType)
    """

    for node, pos in process_next_generator(expression, position):
        # node = equation[pos]
        if isinstance(node, Plus) and is_simple_plus(node):
            return (pos, OperationType.plus)
        elif isinstance(node, Times) and is_simple_times(node):
            return (pos, OperationType.times)

    # [] = position of root
    return ([], OperationType.none)

def process_next_generator(expr, position=[]):
    """ Yields subexpressions to process next.

    Yields (position, expr). Positions are returned in post-order.
    (see https://en.wikipedia.org/wiki/Tree_traversal)

    This funciton is almost the same as
    BaseExpression.iterate_postorder_with_positions(). The difference is that
    this function does not visit the operands of simple sums
    (is_simple_plus(expr) is True). The reason is that those sums can be
    processed directly, so there is no need to yield subexpressions.
    iterate_postorder_with_positions would yield products of the form
    scalar * matrix that are acceptable in simple products.

    # TODO the same could be used for symmetric products
    """
    if is_simple_plus(expr):
        yield (expr, position)

    if not isinstance(expr, Symbol):
        for n, operand in enumerate(expr.operands):
            position_copy = copy.deepcopy(position)
            position_copy.append(n)
            yield from process_next_generator(operand, position_copy)
    yield (expr, position)

def inverse_type(expr):
    """Infers the ExpressionType of the inverse expr."""

    first_operand = expr.operands[0]

    op_type = None
    if isinstance(first_operand, Times):
        op_type = OperationType.times
    elif isinstance(first_operand, Plus):
        op_type = OperationType.plus

    if isinstance(first_operand, Symbol):
        if not first_operand.has_property(properties.ADMITS_FACTORIZATION):
            # if expr is an inverse of a symbol but has a property
            # that makes factorization unnecessary, skip this inverse
            #print(" ".join(["don't factor:", str(first_operand)]) )
            return (ExpressionType.simple_inverse_no_factor, op_type)
        else:
            # otherwise, use factorization
            #print(" ".join(["factor:", str(first_operand)]) )
            return (ExpressionType.simple_inverse, op_type)
    else:
        # if expr is an inverse of a compound expression
        for tmp_expr, pos in expr.preorder_iter():
            # search for operands that could be factored
            if isinstance(tmp_expr, Matrix):
                if first_operand.has_property(properties.ADMITS_FACTORIZATION):
                    # if there is one, use factorizations
                    #print(" ".join(["compound, factor because of:", str(tmp_expr)]) )
                    return (ExpressionType.compound_inverse, op_type)
                else:
                    # if there are none, use reductions
                    continue
                   
        # When we end up here, we know that there are no operands to factor in
        # this compund inverse.
        return (ExpressionType.compound_inverse_no_factor, op_type)

def identify_inverses_eqns(equations):
    """Identifies an innermost inverse in equations."""
    use_factorizations = False
    use_reductions = False

    for eqn_idx, equation in enumerate(equations):
        # always starting at right-hand side
        # inverse_positions uses postorder, so we look at deeper inverses first
        for inv_pos in inverse_positions(equation.rhs, [1]):
            node = equation[inv_pos]
            type = inverse_type(node)
            if type == ExpressionType.simple_inverse_no_factor:
                continue
            elif type == ExpressionType.simple_inverse:
                use_factorizations = True
                return (eqn_idx, inv_pos, use_factorizations, use_reductions)
            elif type == ExpressionType.compound_inverse_no_factor:
                use_reductions = True
                return (eqn_idx, inv_pos, use_factorizations, use_reductions)
            elif type == ExpressionType.compound_inverse:
                use_factorizations = True
                use_reductions = True
                return (eqn_idx, inv_pos, use_factorizations, use_reductions)

    return (None, None, False, False)


def identify_times(equations):
    for n, equation in enumerate(equations):

        eqn_idx = n
        for pos in type_positions(Times, equation.rhs, [1]):
            node = equation[pos]
            if is_simple_times(node):
                return (eqn_idx, pos)

    return (None, None)

def identify_plus(equations):
    for n, equation in enumerate(equations):

        eqn_idx = n
        for pos in type_positions(Plus, equation.rhs, [1]):
            node = equation[pos]
            if is_simple_plus(node):
                return (eqn_idx, pos)
                print(equation)

    return (None, None)    

def is_scalar(expr):
    """Tests if expr only consists of scalars.

    Note:
        This function may return False even expr is a scalar. Because of
        products with vectors, an expression can be scalar even though it does
        not exclusively consist of scalars.
    """
    return all(e.has_property(properties.SCALAR) for e in expr.iterate_preorder())

def is_simple_plus(expr):
    if isinstance(expr, Symbol):
        return False
    return len(expr.operands)>1 and all(is_simple_summand(operand) for operand in expr.operands)

def is_simple_summand(expr):
    if _is_simple_summand(expr):
        return True
    elif isinstance(expr, Times) and len(expr.operands) == 2 and any(isinstance(operand, Scalar) for operand in expr.operands):
        return True
    else:
        return False

def _is_simple_summand(expr):
    if isinstance(expr, Symbol):
        return True
    elif isinstance(expr, Transpose) and isinstance(expr.operand, Symbol):
        return True
    elif isinstance(expr, ConjugateTranspose) and isinstance(expr.operand, Symbol):
        return True
    else:
        return False

def is_simple_times(expr):
    """Test if expr is sufficiently simple for matrix chain rule."""
    return len(expr.operands)>1 and all([_is_simple_times(operand) for operand in expr.operands])

def _is_simple_times(expr):
    # TODO what is missing?
    if isinstance(expr, Symbol):
        return True
    elif isinstance(expr, Transpose) and isinstance(expr.operand, Symbol):
        return True
    elif isinstance(expr, ConjugateTranspose) and isinstance(expr.operand, Symbol):
        return True
    elif isinstance(expr, Inverse) and isinstance(expr.operand, Symbol) and not expr.operand.has_property(properties.ADMITS_FACTORIZATION):
        return True
    elif isinstance(expr, InverseTranspose) and isinstance(expr.operand, Symbol) and not expr.operand.has_property(properties.ADMITS_FACTORIZATION):
        return True
    else:
        return False

def inverse_positions(expr, position=[]):
    """Returns positions of all inverses (including InverseTranspose, …).

    Positions are return in postorder
    (see https://en.wikipedia.org/wiki/Tree_traversal)
    """
    if isinstance(expr, Symbol):
        return

    for n, operand in enumerate(expr.operands):
        position_copy = copy.deepcopy(position)
        position_copy.append(n)
        yield from inverse_positions(operand, position_copy)

    if isinstance(expr, Inverse) or isinstance(expr, InverseTranspose) or isinstance(expr, InverseConjugate) or isinstance(expr, InverseConjugateTranspose):
        yield position
