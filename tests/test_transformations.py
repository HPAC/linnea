from hypothesis import given

from linnea.algebra.transformations import simplify, transpose, invert, invert_transpose

from .generate import expression_strategy


@given(expression_strategy)
def test_transpose_is_self_inverse(expr):
    expr = simplify(expr)
    assert transpose(transpose(expr)) == expr

@given(expression_strategy)
def test_invert_is_self_inverse(expr):
    expr = simplify(expr)
    assert invert(invert(expr)) == expr

@given(expression_strategy)
def test_invert_transpose_is_self_inverse(expr):
    expr = simplify(expr)
    assert invert_transpose(invert_transpose(expr)) == expr

@given(expression_strategy)
def test_transpose_and_invert_order_does_not_matter(expr):
    expr = simplify(expr)
    assert invert(transpose(expr)) == transpose(invert(expr))

@given(expression_strategy)
def test_invert_and_transpose_same_as_invert_transpose(expr):
    expr = simplify(expr)
    assert invert(transpose(expr)) == invert_transpose(expr)

@given(expression_strategy)
def test_transpose_and_invert_same_as_invert_transpose(expr):
    expr = simplify(expr)
    assert transpose(invert(expr)) == invert_transpose(expr)


