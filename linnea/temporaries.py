from .algebra import expression as ae
from .algebra.properties import Property as properties
from .algebra import transformations as at
from .algebra import representations as ar

import copy
import itertools

import matchpy

_counter = 0

def get_identifier():
    """Returns a unique string to indentify temporaries."""
    global _counter
    _counter += 1
    return str(_counter)

# _table_of_temporaries stores temporaries based on the expressions
# they are equivalent to
# Example:
# {Times(A, B, C): tmp1}
_table_of_temporaries = dict()

# For each temporary, it is stored to which expression it is equivalent, (for
# example tmp1 = A + B). This information is then used to reuse temporaries.
_equivalent_expressions = dict()

# self.table_of_factors stores matches of factorization kernels and
# their corresponding output operands. Thus, the same output
# operands can be used for the same input operands.
# This allows to easily remove redundant leaves.
#
# self.table_of_factors is a dictionary of dictionaries:
# keys for self.table_of_factors are kernel IDs.
#   for each ID, there is a dict mapping the matched expression to operands.
# Example:
# {42:
#   {S: L1, P: L2}
#  }
_table_of_factors = dict()

def clear():
    global _table_of_temporaries, _equivalent_expressions, _table_of_factors, _counter
    _table_of_temporaries.clear()
    _equivalent_expressions.clear()
    _table_of_factors.clear()
    _counter = 0

# @profile
def create_tmp(expr, set_equivalent, equiv_expr=None, _properties=None):
    # The value of set_equivalent does not matter if equiv_expr is not None.


    if equiv_expr:
        # This allows to explicitly specify the equivalent expression for expr.
        # If there are temporaries, they are replaced (as a rule of thumb, this
        # should always be done if there are any)
        set_equivalent = True
        # print(equiv_expr)
        equiv_expr = _get_equivalent(equiv_expr)
        # print(equiv_expr)
    elif set_equivalent:
        # print(expr)
        equiv_expr = _get_equivalent(expr)

    # REMEMBER: If it ever turns out to be a problem to use strings here,
    # store real expressions in a sorted list and use bisect to find them.

    # print(expr, set_equivalent, equiv_expr, _properties)

    try:
        tmp = _table_of_temporaries[equiv_expr]
    except KeyError:
        name = "tmp{}".format(get_identifier())
        size = expr.size

        # rows, columns = size
        # if rows != 1 and columns != 1:
        #     tmp = ae.Matrix(name, size)
        # elif rows == 1 and columns == 1:
        #     tmp = ae.Scalar(name)
        # else:
        #     tmp = ae.Vector(name, size)
        if expr.has_property(properties.SCALAR):
            tmp = ae.Scalar(name)
        elif expr.has_property(properties.VECTOR):
            tmp = ae.Vector(name, size)
        elif expr.has_property(properties.MATRIX):
            tmp = ae.Matrix(name, size)

        tmp.indices = expr.indices
        tmp.bandwidth = expr.bandwidth
        if isinstance(expr, ae.Operator) and expr.arity == matchpy.Arity.unary:
            tmp.factorization_labels = expr.factorization_labels

        if set_equivalent:
            # This is effectively a two-way dictionary (or bijective mapping)
            _equivalent_expressions[name] = equiv_expr
            _table_of_temporaries[equiv_expr] = tmp

    # for special properties
    if _properties:
        for prop in _properties:
            tmp.set_property(prop)

    return tmp


# @profile
def _get_equivalent(expr):
    # TODO: This is horrible. Is there any way to fix it?
    new_expr = _flatten_equivalent(expr)
    return ar.to_SOP(at.simplify(ar.to_SOP(new_expr))) if new_expr is not expr else expr

def _flatten_equivalent(expr):
    if isinstance(expr, ae.Symbol):
        return _equivalent_expressions.get(expr.name, expr)
    elif isinstance(expr, ae.Operator):
        new_operands = []
        changed = False
        for operand in expr.operands:
            new_operand = _flatten_equivalent(operand)
            if new_operand is not operand:
                changed = True
            new_operands.append(new_operand)
        return type(expr)(*new_operands) if changed else expr
    return expr



def set_equivalent(expr_before, expr_after):
    """Modifies the mapping of equivalent expressions and temporaries.

    When using matrix factorizations, it is not possible to use equivalent
    expression. To be able to merge branches and avoid redundant computations,
    this function does the following. For the two input expressions, it computes
    the diff paths. For the subexpression of expr_before, a temporary is
    created. For the subexpression of expr_after, its _table_of_temporaries
    entry is set to that temporary. This has the following effect: When an
    expression A*B is computed, eventually, the tables will contain the mapping
    A*B <-> tmp1. Let us now assume the QR factorizatio is applied to B; we
    obtain A*Q*R. Normally, when this expression is computed, another temporary
    (not tmp1) is generated for this expression. By calling
    set_equivalent(A*B, A*Q*R), an additional entry A*Q*R -> tmp1 is added to
    the tables. Thus, once A*Q*R is fully computed, the temporary tmp1 will be
    used for it and the corresponding branches can be merged.

    TL;DR: Adds the following entries to the tables.
    expr_before <-> tmp
    expr_after -> tmp

    Args:
        expr_before (Expression)
        expr_after (Expression)
    """

    # return
    path_before, path_after = expr_diff(expr_before, expr_after)
    # subexpr_before =
    # subexpr_after =
    tmp = create_tmp(expr_before[path_before], True)
    _table_of_temporaries[_get_equivalent(expr_after[path_after])] = tmp

    set_equivalent_upwards(expr_before, expr_after)

def set_equivalent_upwards(expr_before, expr_after):
    """Solves a problem too complicated to be described in one line.

    Under certain circumstances, using set_equivalent causes problems. Let us
    assume we have an expression expr = expr1*expr2. expr1 is computed in many
    different ways, using factorizations. Thus, there will be plenty of entries
    for expr1 in the tables, but as soon as expr1 is fully computed, it will be
    replaced with the same temporary. If a factorization was also applied to
    expr2, the problem is the following: _get_equivalent will replace expr1 with
    the expression that it originally was, without any factorizations applied.
    The resulting expr may not be in the tables, because it is a mixture of
    subexpressions with and without facorizations for which there is no mapping
    to the 'canonical' temporary yet, so a new temporary is generated. This
    functions makes sure that this does not happen by doing the following: After
    transformations that may cause this problem (probably all of them), we
    compute the diff paths of expr_before and expr_after (the root expressions).
    Now, we look for the subexpression of expr_before with the longest prefix of
    its diff path that is in the tables (it should always exist, in the worst
    case, it is the root expression). Then, we take the prefix of the same
    length for expr_after and add entries for this expression to the tables.

    Reminder: The described problem happens in example 63 when removing this
    function from TR_matrix_chain.

    Args:
        expr_before (Expression)
        expr_after (Expression)
    """

    # return
    path_before, path_after = expr_diff(expr_before, expr_after)
    # print("here", path_before)
    for i in reversed(range(len(path_before))):
        # print(path_before[:i])
        equiv_subexpr_before = _get_equivalent(expr_before[path_before[:i]])
        if equiv_subexpr_before in _table_of_temporaries:
            equiv_subexpr_after = _get_equivalent(expr_after[path_after[:i]])
            if equiv_subexpr_after not in _table_of_temporaries:
                # print(str(_get_equivalent(expr_before[path_before[:i]])))
                # print(str(_get_equivalent(expr_after[path_after[:i]])))
                tmp = create_tmp(equiv_subexpr_before, True)
                _table_of_temporaries[equiv_subexpr_after] = tmp
                break


def expr_diff(e1, e2, p1=[], p2=[]):
    """Computes the diff paths of two expressions.

    For two expressions e1!=e2, there exist two path p1, p2 such that by
    replacing the subexpressions e1[p1] and e2[p2] with the same symbol, one
    optains two identical expressions. To put it differently, p1 and p2 identify
    those positions in e1 and e2, respectively, where those two expressions
    differ. In the "worst" case where e1 and e2 are completely different
    expressions, p1 and p2 will be the empty paths. For syntactic expressions,
    it always holds that p1==p2.

    Args:
        e1 (Expression)
        e2 (Expression)
        p1 (list, optional): By default, this is the empty path. Only relevant
            for recursive calls.
        p2 (list, optional): By default, this is the empty path. Only relevant
            for recursive calls.
    """
    if e1 == e2:
        return None

    if type(e1) == type(e2) and isinstance(e1, ae.Operator) and len(e1.operands) == len(e2.operands):
        if e1.commutative:
            """
            - go over all permutions of the operands of one expr
            - look at all permuations
              - if there is one with exactly one mismatch, go deeper (we
                have found what we are looking for)
              - if there is one with zero mismatches, abort (should never happen because of e1!=e2)
              - if all have more than one, we have found the deepest node
            """
            for permuted_operands in itertools.permutations(enumerate(e2.operands)):
                # print(permuted_operands)
                rets = []
                for (pos1, op1), (pos2, op2) in zip(enumerate(e1.operands), permuted_operands):
                    ret = expr_diff(op1, op2, p1+[pos1], p2+[pos2])
                    # print(op1, op2, ret)
                    if ret:
                        rets.append(ret)
                # print(rets)
                if len(rets)==1:
                    return rets[0]
        else:
            # at most one pair can returns something, otherwise, the current
            # positions are what we are looking for
            rets = []
            for pos, (op1, op2) in enumerate(zip(e1.operands, e2.operands)):
                ret = expr_diff(op1, op2, p1+[pos], p2+[pos])
                # print(op1, op2, ret)
                if ret:
                    rets.append(ret)
            if len(rets)==1:
                return rets[0]
    return (p1, p2)