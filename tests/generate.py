from random import random
import math

from hypothesis import assume, given, find
import hypothesis.strategies as st

import linnea.algebra.expression as ae
from linnea.algebra.equations import Equations
from linnea.algebra.transformations import simplify, transpose

__all__ = ['expression_strategy']


def func_wrap_strategy(args, func):
    min_size = func.arity[0] or 1
    max_size = func.arity[0] if func.arity[1] else 4
    return st.lists(
        args, min_size=min_size, max_size=max_size).map(lambda a: func(*a))


def expression_recurse_strategy(args):
    return (func_wrap_strategy(args, ae.Times) |
            func_wrap_strategy(args, ae.Plus) |
            func_wrap_strategy(args, ae.Transpose) |
            func_wrap_strategy(args, ae.Inverse) |
            func_wrap_strategy(args, ae.InverseTranspose))


class NeedQuadratic(Exception):
    pass


def random_dimension():
    # random number biased towards lower values
    return 1 + math.floor(random()**2 * 20)


def randomize_sizes(expr, rows=None, cols=None):
    can_change = False
    if rows is None:
        rows = random_dimension()
    if cols is None:
        can_change = True
        cols = random_dimension()
    if isinstance(expr, ae.Symbol):
        if rows == 1 or cols == 1:
            if rows == cols:
                return ae.Scalar(expr.name.replace('M', 'a'))
            else:
                return ae.Vector(expr.name.replace('M', 'v'), (rows, cols))
        return ae.Matrix(expr.name, (rows, cols))
    try:
        if isinstance(expr, ae.Times):
            operands = []
            new_rows = rows
            for operand in expr.operands[:-1]:
                new_cols = random_dimension()
                try:
                    operands.append(
                        randomize_sizes(operand, new_rows, new_cols))
                    new_rows = new_cols
                except NeedQuadratic:
                    operands.append(
                        randomize_sizes(operand, new_rows, new_rows))
            try:
                operands.append(
                    randomize_sizes(expr.operands[-1], new_rows, cols))
            except NeedQuadratic:
                operands.append(
                    randomize_sizes(expr.operands[-1], new_rows, new_rows))
            return ae.Times(*operands)
        if isinstance(expr, ae.Transpose):
            return ae.Transpose(randomize_sizes(expr.operand, cols, rows))
        if isinstance(expr, (ae.Inverse, ae.InverseTranspose)):
            if rows != cols:  # Inversion needs a quadratic matrix or a scalar
                raise NeedQuadratic()
            return type(expr)(randomize_sizes(expr.operand, rows, rows))
        if isinstance(expr, ae.Operator):
            return type(expr)(*(randomize_sizes(operand, rows, cols)
                                for operand in expr.operands))
        assert False, "Unreachable"
    except NeedQuadratic:
        if not can_change:
            raise
        return randomize_sizes(expr, rows, rows)


_n = 0


def next_matrix():
    global _n
    _n += 1
    return ae.Matrix('M{}'.format(_n), size=(5, 5))


expression_base_strategy = st.builds(lambda _: next_matrix(), st.uuids())
expression_strategy = st.tuples(
    st.recursive(
        expression_base_strategy, expression_recurse_strategy, max_leaves=20),
    st.
    random_module()  # needed to seed the random module for randomize_sizes()
).map(lambda x: randomize_sizes(x[0]))
