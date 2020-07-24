from ..algebra.expression import Equal, Matrix, Vector, Scalar
from ..algebra.properties import Property
from ..algebra.equations import Equations

from .. import temporaries

from .graph import properties as gp

_counter = 0


def add_expression(expr, properties):
    # TODO what happens if the same expression is added a second time
    # - with the same properties?
    # - with different properties?

    global _counter
    name = "".join(["SP", str(_counter)])
    _counter +=1
    lhs = None
    if expr.has_property(Property.MATRIX):
        lhs = Matrix(name, expr.size)
    elif expr.has_property(Property.VECTOR):
        lhs = Vector(name, expr.size)
    elif expr.has_property(Property.SCALAR):
        lhs = Scalar(name)

    graph = gp.PropertyGraph(Equations(Equal(lhs, expr)))
    graph.derivation()

    # graph.write_output(code=False,
    #                    derivation=False,
    #                    output_name=name,
    #                    experiment_code=False,
    #                    k_best=False,
    #                    algorithms_limit=100,
    #                    pruning_factor=1.0,
    #                    graph=True,
    #                    # graph_style=linnea.config.GraphStyle.full,
    #                    no_duplicates=True)

    for node in graph.nodes:
        expr = node.equations[-1].rhs
        tmp = temporaries.create_tmp(expr)
        for prop in properties:
            tmp.set_property(prop)
