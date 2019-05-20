from ..algebra.expression import Equal, Matrix, Vector, Scalar
from ..algebra.properties import Property as properties
from ..algebra.equations import Equations

from .. import temporaries

from .graph import properties as gp

_counter = 0


def add_expression(expr, _properties):
    # TODO what happens if the same expression is added a second time
    # - with the same properties?
    # - with different properties?

    global _counter
    name = "".join(["SP", str(_counter)])
    _counter +=1
    lhs = None
    if expr.has_property(properties.MATRIX):
        lhs = Matrix(name, expr.size)
    elif expr.has_property(properties.VECTOR):
        lhs = Vector(name, expr.size)
    elif expr.has_property(properties.SCALAR):
        lhs = Scalar(name)

    graph = gp.PropertyGraph(Equations(Equal(lhs, expr)))
    graph.derivation()

    for node in graph.nodes:
        expr = node.equations[-1].rhs
        tmp = temporaries.create_tmp(expr, True, _properties=_properties)
