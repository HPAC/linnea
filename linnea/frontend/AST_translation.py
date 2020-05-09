from tatsu.walkers import NodeWalker

from ..algebra import expression as ae
from ..algebra.equations import Equations
from ..algebra.properties import Property


class LinneaWalker(NodeWalker):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._variables = dict()
        self._symbols = dict()
        self._equations = []
        self.symbolic_operand_sizes = dict()

    @property
    def equations(self):
        return Equations(*self._equations)

    def walk_object(self, node):
        assert False, "This should never be reached unless the grammar is changed."

    def walk_Model(self, node):
        for var in node.vars:
            self.walk(var)
        for symbol in node.symbols:
            self.walk(symbol)
        for equation in node.equations:
            self.walk(equation)

    def walk_Size(self, node):
        self._variables[node.name] = int(node.value)

    def _set_symbol(self, symbol, node):
        for prop in node.properties:
            symbol.set_property(Property(prop))
        self._symbols[node.name] = symbol

    def walk_Matrix(self, node):
        size = (self._variables[node.dims.rows], self._variables[node.dims.columns])
        self.symbolic_operand_sizes[node.name] = (node.dims.rows, node.dims.columns)
        self._set_symbol(ae.Matrix(node.name, size), node)

    def walk_RowVector(self, node):
        length = self._variables[node.dims.length]
        self.symbolic_operand_sizes[node.name] = (node.dims.length)
        self._set_symbol(ae.Vector(node.name, (1, length)), node)

    def walk_ColumnVector(self, node):
        length = self._variables[node.dims.length]
        self.symbolic_operand_sizes[node.name] = (node.dims.length)
        self._set_symbol(ae.Vector(node.name, (length, 1)), node)

    def walk_IdentityMatrix(self, node):
        size = (self._variables[node.dims.rows], self._variables[node.dims.columns])
        self.symbolic_operand_sizes[node.name] = (node.dims.rows, node.dims.columns)
        self._symbols[node.name] = ae.IdentityMatrix(*size)

    def walk_ZeroMatrix(self, node):
        size = (self._variables[node.dims.rows], self._variables[node.dims.columns])
        self.symbolic_operand_sizes[node.name] = (node.dims.rows, node.dims.columns)
        self._symbols[node.name] = ae.ZeroMatrix(*size)

    def walk_Scalar(self, node):
        self._set_symbol(ae.Scalar(node.name), node)

    def walk_Number(self, node):
        return ae.ConstantScalar(float(node.value))

    def walk_Plus(self, node):
        return ae.Plus(self.walk(node.left), self.walk(node.right))

    def walk_Times(self, node):
        return ae.Times(self.walk(node.left), self.walk(node.right))

    def walk_Subtract(self, node):
        return ae.Plus(self.walk(node.left), ae.Times(ae.ConstantScalar(-1.0), self.walk(node.right)))

    def walk_Transpose(self, node):
        return ae.Transpose(self.walk(node.arg))

    def walk_Inverse(self, node):
        return ae.Inverse(self.walk(node.arg))

    def walk_Minus(self, node):
        return ae.Times(ae.ConstantScalar(-1.0), self.walk(node.arg))

    def walk_Symbol(self, node):
        return self._symbols[node.name]

    def walk_Equation(self, node):
        self._equations.append(ae.Equal(self.walk(node.lhs), self.walk(node.rhs)))

# file_name = "input_test.la"
# with open(file_name, "r") as input_file:
#     my_parser = parser.LinneaParser()

#     ast = my_parser.parse(input_file.read(), rule_name = "model")

#     print(json.dumps(grako.util.asjson(ast.assignments), indent=2))

#     ast_translator = ASTTranslator(ast)
#     print(ast_translator.equations)
