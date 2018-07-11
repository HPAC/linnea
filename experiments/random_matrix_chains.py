
import random
import math
import numpy


import itertools
from linnea.algebra import expression as ae
from linnea.algebra.equations import Equations
from linnea.algebra.transformations import simplify
from linnea.algebra.properties import Property as properties
from linnea.frontend.export import export
from linnea.code_generation import experiments as cge
from linnea.utils import window


def matrix_size():
    return random.randrange(50, 2001, 50)

_n = 0

def generate_operand(rows, columns):
    global _n
    _n += 1
    if columns == rows:
        if rows == 1:
            return # scalar

        operand = ae.Matrix('M{}'.format(_n), size=(rows, columns))
        operand.set_property(properties.FULL_RANK)
        # include "no property"
        if random.random() > 0.25:
            operand.set_property(random.choices([properties.DIAGONAL, properties.LOWER_TRIANGULAR, properties.UPPER_TRIANGULAR, properties.SYMMETRIC, properties.SPD])[0])
        return operand
    elif columns == 1:
        return ae.Vector('v{}'.format(_n), size=(rows, columns))
    elif rows == 1:
        return ae.Vector('v{}'.format(_n), size=(rows, columns))
    else:
        operand = ae.Matrix('M{}'.format(_n), size=(rows, columns))
        operand.set_property(properties.FULL_RANK)
        return operand

def matrix_chain_generator():
    while True:
        # for _ in range(20):
        #     expression_strategy.example()
        # generate_matrix_chain()

        length = random.randrange(4, 10)
        sizes = []

        if  0.95 <= random.random():
            sizes.append(1)
        else:
            sizes.append(matrix_size())

        for i in range(length):
            rand = random.random()
            if  0.6 <= rand < 0.95:
                # square
                sizes.append(sizes[i])
            elif 0.95 <= rand:
                # vector
                sizes.append(1)
            else:
                # non-square
                sizes.append(matrix_size())

        # print(sizes)
        operands = []
        restart = False
        for rows, columns in window(sizes):
            if rows == columns:
                if rows == 1:
                    restart = True
                    break
                    # operands.append(None) # error
                op = random.choices([ae.Identity, ae.Transpose, ae.Inverse, ae.InverseTranspose], weights=[2, 1, 1, 1])[0]
                operands.append(op(generate_operand(rows, columns)))
            elif rows == 1:
                operands.append(ae.Transpose(generate_operand(columns, rows)))
            elif columns == 1:
                operands.append(generate_operand(rows, columns))
            else:
                op = random.choices([ae.Identity, ae.Transpose], weights=[3, 1])[0]
                if op == ae.Transpose:
                    operands.append(op(generate_operand(columns, rows)))
                else:
                    operands.append(op(generate_operand(rows, columns)))
        if restart:
            continue

        # print(operands)
        expr = ae.Times(*operands)
        expr = simplify(expr)
        # print(expr)

        # print(expr.size)
        # print((sizes[0], sizes[-1]))
        if expr.has_property(properties.MATRIX):
            lhs = ae.Matrix('X'.format(_n), size=(sizes[0], sizes[-1]))
        elif expr.has_property(properties.VECTOR):
            lhs = ae.Vector('x'.format(_n), size=(sizes[0], sizes[-1]))

        yield Equations(ae.Equal(lhs, expr))

