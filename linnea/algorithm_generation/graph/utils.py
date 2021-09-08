from ...algebra.expression import Symbol, Matrix, Scalar, Inverse, \
                                  Times, Plus, Operator, \
                                  Transpose, ConjugateTranspose, \
                                  InverseConjugate, \
                                  InverseConjugateTranspose

from ...algebra.transformations import transpose, invert, undistribute_inverse, \
                                       admits_undistribution
from ...algebra.representations import to_POS
from ...algebra import equations as aeq
from ...algebra.properties import Property

from ...utils import is_inverse, is_transpose

from ..utils import is_blocked

from enum import Enum, unique

import copy
import itertools
import operator
import heapq

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
    unary = 3
    none = 4


@unique
class GenerationStep(Enum):
    factorizations = 1
    kernels = 2
    tricks = 3
    CSE = 4
    prune = 5
    merge = 6


def find_operands_to_factor(equations, eqn_idx=None):
    """Finds all operands to factor.

    Finds all operands that may require factorizations. Those are operands that
    have the property ADMITS_FACTORIZATION and appear within an inverse.

    Args:
        equations (Equations): The equations that are searched.

    Returns:
        set: A set of operands.
    """
    # this is independent of variants (most likely)

    eqn_indices = None
    if eqn_idx is not None:
        eqn_indices = [eqn_idx]
    else:
        eqn_indices = range(len(equations))

    operands_to_factor = set()
    for _eqn_idx in eqn_indices:
        equation = equations[_eqn_idx]
        for inv_expr, inv_pos in inverse_positions(equation.rhs, (1,)):
            for expr, pos in inv_expr.preorder_iter():
                if isinstance(expr, Symbol) and expr.has_property(Property.ADMITS_FACTORIZATION):
                    operands_to_factor.add(expr)

    """Operands that have not been computed yet can not be factored. An operand
    tmp is not computed yet if there is an equation tmp = expr where expr is not
    a Symbol.
    """
    for _eqn_idx in eqn_indices:
        equation = equations[_eqn_idx]
        if equation.lhs in operands_to_factor and not isinstance(equation.rhs, Symbol):
            operands_to_factor.remove(equation.lhs)

    return operands_to_factor


def find_explicit_symbol_inverse(expr, position=(), predecessor=None):

    if is_inverse(expr) and isinstance(expr.operand, Symbol) and not isinstance(predecessor, Times):
        yield (expr, position)

    if isinstance(expr, Operator):
        for n, operand in enumerate(expr.operands):
            new_position = position + (n,)
            yield from find_explicit_symbol_inverse(operand, new_position, expr)


# @profile
def generate_representations(equations, eqn_idx=None):
    """Generates different representations of equations.

    This function yields the following representations of equations (if they are
    unique), in the order given below:
        1. POS variant (left first).
        2. POS variant (right first).
        3. UI variant.
        4. Equations as passed to the function.

    If an eqn_idx is given, different representations are only generated for the
    respective equation.

    Args:
        equations (linnea.algebra.equations.Equations): Input equations.

        eqn_idx (int, optional): Index of equation to generate variants for.

    Yields:
        linnea.algebra.equations.Equations: Different representations of
            equations.

    """

    # TODO what about combinations of POS and undistribute inverse?
    # try all combinations? yes, but only if the same equations are concerned

    yielded_variants = set()

    if eqn_idx is None:
        eqn_indices = range(len(equations))
    else:
        eqn_indices = [eqn_idx]

    # TODO combine this with the other loop
    POS_candidates = set()
    for _eqn_idx in eqn_indices:
        equation = equations[_eqn_idx]
        for node, pos in equation.preorder_iter():
            if isinstance(node, Plus):
                for operand in node.operands:
                    if isinstance(operand, Times):
                        POS_candidates.add(_eqn_idx)
                        break

    if POS_candidates:
        new_equations_left = []
        for _eqn_idx, equation in enumerate(equations):
            new_equation_left = equation
            if _eqn_idx in POS_candidates:
                new_equation_left = to_POS(new_equation_left, "l")
            new_equations_left.append(new_equation_left)

        temp_eqns = aeq.Equations(*new_equations_left)
        if temp_eqns not in yielded_variants:
            yielded_variants.add(temp_eqns)
            yield temp_eqns

        new_equations_right = []
        for _eqn_idx, equation in enumerate(equations):
            new_equation_right = equation
            if _eqn_idx in POS_candidates:
                new_equation_right = to_POS(new_equation_right, "r")
            new_equations_right.append(new_equation_right)

        temp_eqns = aeq.Equations(*new_equations_right)
        if temp_eqns not in yielded_variants:
            yielded_variants.add(temp_eqns)
            yield temp_eqns


    undist_inv_candidates = set()
    for _eqn_idx in eqn_indices:
        equation = equations[_eqn_idx]
        for node, pos in equation.preorder_iter():
            if isinstance(node, Times):
                number_of_inverses = 0
                for operand in node.operands:
                    if admits_undistribution(operand):
                        number_of_inverses += 1
                    if number_of_inverses > 1:
                        undist_inv_candidates.add(_eqn_idx)
                        break

    if undist_inv_candidates:
        new_equations = []
        for _eqn_idx, equation in enumerate(equations):
            new_equation = equation
            if _eqn_idx in undist_inv_candidates:
                new_equation = undistribute_inverse(new_equation)
            new_equations.append(new_equation)

        temp_eqns = aeq.Equations(*new_equations)
        if temp_eqns not in yielded_variants:
            yielded_variants.add(temp_eqns)
            yield temp_eqns

    if equations not in yielded_variants:
        yield equations



def process_next(equation):
    """Finds a subexpression to process next in the generation.

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
    for expr, pos in inverse_positions(equation.rhs, (1,)):
        expr_type, op_type = inverse_type(expr)
        if expr_type == ExpressionType.simple_inverse_no_factor:
            continue
        else:
            return (pos, expr_type, op_type)


    for expr, pos in process_next_generator(equation.rhs, (1,)):
        # expr = equation[pos]
        # print("here", expr, expr.size)
        # print(list((n, n.size) for n in expr.iterate_preorder()))
        # if is_scalar(expr):
            # print("here")
        if isinstance(expr, Plus) and is_simple_plus(expr):
            return (pos, ExpressionType.simple, OperationType.plus)
        elif isinstance(expr, Times) and is_simple_times(expr):
            return (pos, ExpressionType.simple, OperationType.times)

    # (1,) = position of right-hand side
    return ((1,), ExpressionType.none, OperationType.none)


def process_next_simple(expression):
    """Finds a subexpression to process next in the generation.

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

    if (is_inverse(expression.rhs) or is_transpose(expression.rhs)) and isinstance(expression.rhs.operand, Symbol):
        return ((1,), OperationType.unary)    

    for node, pos in process_next_generator(expression):
        # node = equation[pos]
        if isinstance(node, Plus) and is_simple_plus(node):
            return (pos, OperationType.plus)
        elif isinstance(node, Times) and is_simple_times(node):
            return (pos, OperationType.times)

    # () = position of root
    return ((), OperationType.none)


def process_next_generator(expr, position=()):
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
            new_position = position + (n,)
            yield from process_next_generator(operand, new_position)
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
        if not first_operand.has_property(Property.ADMITS_FACTORIZATION):
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
                if first_operand.has_property(Property.ADMITS_FACTORIZATION):
                    # if there is one, use factorizations
                    #print(" ".join(["compound, factor because of:", str(tmp_expr)]) )
                    return (ExpressionType.compound_inverse, op_type)
                else:
                    # if there are none, use other kernels
                    continue
                   
        # When we end up here, we know that there are no operands to factor in
        # this compund inverse.
        return (ExpressionType.compound_inverse_no_factor, op_type)


def identify_inverses_eqns(equations):
    """Identifies an innermost inverse in equations."""
    use_factorizations = False
    use_kernels = False

    for eqn_idx, equation in enumerate(equations):
        # always starting at right-hand side
        # inverse_positions uses postorder, so we look at deeper inverses first
        for inv_expr, inv_pos in inverse_positions(equation.rhs, (1,)):
            type = inverse_type(inv_expr)
            if type == ExpressionType.simple_inverse_no_factor:
                continue
            elif type == ExpressionType.simple_inverse:
                use_factorizations = True
                return (eqn_idx, inv_pos, use_factorizations, use_kernels)
            elif type == ExpressionType.compound_inverse_no_factor:
                use_kernels = True
                return (eqn_idx, inv_pos, use_factorizations, use_kernels)
            elif type == ExpressionType.compound_inverse:
                use_factorizations = True
                use_kernels = True
                return (eqn_idx, inv_pos, use_factorizations, use_kernels)

    return (None, None, False, False)


def is_explicit_inversion(matrix_chain):
    """Tests if matrix_chain is an explicit inversion.

    A matrix chain is an explicit inversion if all operands that are not scalars
    come from the same factorization and at least one operand is inverted. We do
    not require all operands to be inverted because of orthogonal matrices.

    Args:
        matrix_chain (Times): The matrix chain to test.

    Returns:
        bool: True if the chain is an explicit inversion, False otherwise.
    """
    if any(is_inverse(operand) for operand in matrix_chain.operands):
        labels = set()
        for operand in matrix_chain.operands:
            if not operand.has_property(Property.SCALAR):
                if operand.factorization_labels:
                    labels.add(frozenset(operand.factorization_labels))
                else:
                    return False

        return len(labels) == 1
    else:
        return False


def is_dead_end(equations, factored_operands):
    for equation in equations:
        for expr, _ in equation.rhs.preorder_iter():
            if is_inverse(expr) and expr.operand.has_property(Property.ADMITS_FACTORIZATION) and expr.operand in factored_operands:
                return True
            if isinstance(expr, Times) and is_blocked(expr) and not is_explicit_inversion(expr):
                return True
    return False


def is_scalar(expr):
    """Tests if expr only consists of scalars.

    Note:
        This function may return False even expr is a scalar. Because of
        products with vectors, an expression can be scalar even though it does
        not exclusively consist of scalars.
    """
    return all(e.has_property(Property.SCALAR) for e in expr.iterate_preorder())


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
    else:
        # is_transpose cannot be used here because it includes inverse
        return isinstance(expr, (Transpose, ConjugateTranspose)) and isinstance(expr.operand, Symbol)


def is_simple_times(expr):
    """Test if expr is sufficiently simple for matrix chain algorithm."""
    return len(expr.operands)>1 and all([_is_simple_times(operand) for operand in expr.operands])


def _is_simple_times(expr):
    # TODO what is missing?
    if isinstance(expr, Symbol):
        return True
    # order is important here because is_transpose includes invers
    elif is_inverse(expr) and isinstance(expr.operand, Symbol) and not expr.operand.has_property(Property.ADMITS_FACTORIZATION):
        return True
    elif is_transpose(expr) and isinstance(expr.operand, Symbol):
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
        new_position = position + (n,)
        yield from inverse_positions(operand, new_position)

    if is_inverse(expr):
        yield expr, position


class PriorityStack():
    def __init__(self):
        self.heap = []
        self.counter = 0

    def put(self, prio, elem):
        heapq.heappush(self.heap, (prio, -self.counter, elem))
        self.counter += 1

    def get(self):
        prio, _, elem = heapq.heappop(self.heap)
        return prio, elem

    def empty(self):
        return not self.heap
