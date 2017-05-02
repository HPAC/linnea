

# import grako
# import parser
# import sys
# import json

from ..algebra import expression as ae
from ..algebra.equations import Equations
from ..algebra.properties import Property as properties

class ASTTranslator(object):
    """docstring for ASTTranslator"""
    def __init__(self, ast):
        super(ASTTranslator, self).__init__()
        self.ast = ast
        self.variables = dict()
        self.operands = dict()
        self.equations = []
        self.symbolic_operand_sizes = dict()

        self.var_declarations()
        self.op_declarations()
        self.assignments()
        # print(self.equations)
        # print(self.symbolic_operand_sizes)

    def var_declarations(self):
        for declaration in self.ast.var_declarations:
            self.variables[declaration.lhs] = int(declaration.rhs)

    def op_declarations(self):
        for declaration in self.ast.op_declarations:
            if declaration.optype is "Matrix":
                size = (self.variables[declaration.dims.rows], self.variables[declaration.dims.columns])
                self.symbolic_operand_sizes[declaration.name] = (declaration.dims.rows, declaration.dims.columns)
                operand = ae.Matrix(declaration.name, size)
            elif declaration.optype is "RowVector":
                length = self.variables[declaration.dims.length]
                self.symbolic_operand_sizes[declaration.name] = (declaration.dims.length)
                operand = ae.Vector(declaration.name, (1, length))
            elif declaration.optype is "ColumnVector":
                length = self.variables[declaration.dims.length]
                self.symbolic_operand_sizes[declaration.name] = (declaration.dims.length)
                operand = ae.Vector(declaration.name, (length, 1))
            elif declaration.optype is "IdentityMatrix":
                size = (self.variables[declaration.dims.rows], self.variables[declaration.dims.columns])
                self.symbolic_operand_sizes[declaration.name] = (declaration.dims.rows, declaration.dims.columns)
                operand = ae.IdentityMatrix(*size)

            for property in declaration.properties:
                operand.set_property(properties(property))
            self.operands[operand.name] = operand

    def assignments(self):
        equations = []
        for assignment in self.ast.assignments:
            lhs = self.operands[assignment.lhs]
            rhs = self.expression(assignment.rhs)
            equations.append(ae.Equal(lhs, rhs))
        self.equations = Equations(*equations)

    def expression(self, expression):
        if isinstance(expression, str):
            try:
                scalar = float(expression)
            except ValueError:
                return self.operands[expression]
            else:
                return ae.ConstantScalar(scalar)

        args = [self.expression(arg) for arg in expression.args]
        if expression.op == "*":
            return ae.Times(*args)
        elif expression.op == "+":
            return ae.Plus(*args)
        elif expression.op == "-":
            if len(args) == 2:
                return ae.Plus(args[0], ae.Times(ae.ConstantScalar(-1), args[1]))
            elif len(args) == 1:
                return ae.Times(ae.ConstantScalar(-1), args[0])
        elif expression.op == "'" or expression.op == "trans":
            return ae.Transpose(*args)
        elif expression.op == "inv":
            return ae.Inverse(*args)
        else:
            raise NotImplementedError(expression.op)

# file_name = "input_test.la"
# with open(file_name, "r") as input_file:
#     my_parser = parser.LinneaParser()

#     ast = my_parser.parse(input_file.read(), rule_name = "model")

#     print(json.dumps(grako.util.asjson(ast.assignments), indent=2))

#     ast_translator = ASTTranslator(ast)
#     print(ast_translator.equations)
