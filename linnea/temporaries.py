from .algebra import expression as ae
from .algebra.properties import Property
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

equivalence_replacer = matchpy.ManyToOneReplacer()

def clear():
    global _table_of_temporaries, _equivalent_expressions, _table_of_factors, \
           _counter, equivalence_replacer
    _table_of_temporaries.clear()
    _equivalent_expressions.clear()
    _table_of_factors.clear()
    equivalence_replacer = matchpy.ManyToOneReplacer()
    """It is currently not possible to reset _counter here because MatchPy
    requires that variable names are unique. Without unique variable names, if
    the generation is run multiple times in the same program, it
    can happen that the substitute function returns a wrong variable with the
    same name.
    """
    # _counter = 0

def create_tmp(expr):

    tmp = _table_of_temporaries.get(expr)
    if tmp:
        return tmp

    equiv_expr = _flatten_equivalent(expr)
    tmp = _table_of_temporaries.get(equiv_expr)
    if tmp:
        if expr != tmp:
            # Is this necessary? What about the other cases?
            _table_of_temporaries[expr] = tmp
        return tmp

    equiv_expr_normal = ar.to_normalform(equiv_expr)
    tmp = _table_of_temporaries.get(equiv_expr_normal)
    if tmp:
        _table_of_temporaries[equiv_expr] = tmp
        _table_of_temporaries[expr] = tmp
        return tmp

    equiv_expr_normal_replaced = equivalence_replacer.replace(equiv_expr_normal)
    tmp = _table_of_temporaries.get(equiv_expr_normal_replaced)
    if tmp:
        _table_of_temporaries[equiv_expr_normal] = tmp
        _table_of_temporaries[equiv_expr] = tmp
        _table_of_temporaries[expr] = tmp
        return tmp

    # TODO what about another conversion to normal form?

    name = "tmp{}".format(get_identifier())
    size = expr.size

    if expr.has_property(Property.SCALAR):
        tmp = ae.Scalar(name)
    elif expr.has_property(Property.VECTOR):
        tmp = ae.Vector(name, size)
    elif expr.has_property(Property.MATRIX):
        tmp = ae.Matrix(name, size)

    tmp.indices = expr.indices
    tmp.bandwidth = expr.bandwidth
    if isinstance(expr, ae.Operator) and all(operand.factorization_labels for operand in expr.operands):
        # This covers unary operations, but also the case when an operation
        # is not blocked because the result does not admit factorizations.
        tmp.factorization_labels = expr.factorization_labels

    # This is effectively a two-way dictionary (or bijective mapping)
    _equivalent_expressions[name] = equiv_expr_normal_replaced
    _table_of_temporaries[equiv_expr_normal_replaced] = tmp
    _table_of_temporaries[equiv_expr_normal] = tmp
    _table_of_temporaries[equiv_expr] = tmp
    _table_of_temporaries[expr] = tmp

    # print(tmp, expr, equiv_expr)

    return tmp


def get_equivalent(expr):
    return _equivalent_expressions[create_tmp(expr).name]




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
