from ..algebra import expression as ae
from ..algebra.properties import Property as properties


_BINARY = {
    ae.Plus: ' + ',
    ae.Times: ' * ',
}

_UNARY = {
    ae.Transpose: 'trans',
    ae.Inverse: 'inv',
}

_SUPPORTED_PROPERTIES = {
    'Square', 'SPD', 'ColumnPanel', 'RowPanel', 'Diagonal', 'Tridiagonal',
    'Banded', 'LowerTriangular', 'UpperTriangular', 'UnitDiagonal', 'Symmetric',
    'Hessenberg', 'Orthogonal', 'FullRank', 'Non-singular'
}

def export(equations):
    """Export a set of equations to the input string format."""
    # TODO: Add support for names of sizes once they are saved with the symbols
    var_declarations = []
    op_declarations = []
    symbols = set()
    for eq in equations:
        for expr, _ in eq.preorder_iter(lambda e: isinstance(e, ae.Symbol)):
            if expr not in symbols:
                props = ['InOut'] + [p.value for p in expr.properties if p.value in _SUPPORTED_PROPERTIES]
                prop_str = ', '.join(props)
                if isinstance(expr, ae.Matrix):
                    var_declarations.append(
                        '{}_rows = {}'.format(expr.name, expr.size[0]))
                    var_declarations.append(
                        '{}_cols = {}'.format(expr.name, expr.size[1]))
                    if expr.has_property(properties.IDENTITY):
                        var_declarations.append(
                            'IdentityMatrix {0} ({0}_rows, {0}_cols)'.format(expr.name))
                    else:
                        op_declarations.append(
                            'Matrix {0} ({0}_rows, {0}_cols) <{1}>'.format(expr.name, prop_str))
                elif isinstance(expr, ae.Vector):
                    var_declarations.append(
                        '{}_size = {}'.format(expr.name, max(expr.size)))
                    op_declarations.append('{2}Vector {0} ({0}_size) <{1}>'.format(
                        expr.name, prop_str, 'Row' if expr.size[0] == 1 else 'Column'))
                elif isinstance(expr, ae.Scalar):
                    op_declarations.append('Scalar {} <{}>'.format(expr.name, prop_str))
                else:
                    raise TypeError(
                        'Unsupported operand type: {!r}'.format(type(expr)))
                symbols.add(expr)

    assignments = [export_expression(e) for e in equations]

    return '{}\n\n{}\n\n{}'.format(
        '\n'.join(sorted(var_declarations)),
        '\n'.join(sorted(op_declarations)),
        '\n'.join(assignments)
    )


def export_expression(expression):
    """Convert an expression to the input string format."""
    # TODO: Add support for conjugate etc. once they are in the input language
    if isinstance(expression, ae.Symbol):
        return expression.name
    if isinstance(expression, ae.ConstantScalar):
        return str(expression.value)
    if isinstance(expression, ae.Equal):
        return '{} = {}'.format(
            export_expression(expression.lhs),
            export_expression(expression.rhs)
        )
    for op_type, op_str in _BINARY.items():
        if isinstance(expression, op_type):
            return '({})'.format(op_str.join(export_expression(e) for e in expression.operands))
    for op_type, op_str in _UNARY.items():
        if isinstance(expression, op_type):
            return '{}({})'.format(op_str, export_expression(expression.operand))
    if isinstance(expression, ae.InverseTranspose):
        return 'inv(trans({}))'.format(export_expression(expression.operand))
    raise TypeError("Unknown expression type: {!r}".format(type(expression)))
