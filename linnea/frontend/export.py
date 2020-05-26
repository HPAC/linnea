from ..algebra import expression as ae
from ..algebra.properties import Property

from ..utils import dependent_dimensions


_BINARY = {
    ae.Plus: ' + ',
    ae.Times: ' * ',
}

_UNARY = {
    ae.Transpose: 'trans',
    ae.Inverse: 'inv',
}

_SUPPORTED_PROPERTIES = {
    'SPD', 'SPSD', 'Diagonal', 'Permutation', 'Positive',
    'LowerTriangular', 'UpperTriangular', 'UnitDiagonal', 'Symmetric',
    'Orthogonal', 'OrthogonalRows', 'OrthogonalColumns', 'FullRank', 'Non-singular'
}


def export(equations):
    """Export a set of equations to the input string format."""
    dimensions_dict = dict()

    operands = set()
    required_dimensions = set()
    for eq in equations:
        for op, _ in eq.preorder_iter(lambda e: isinstance(e, ae.Symbol) and not isinstance(e, ae.ConstantScalar)):
            operands.add(op)
            if op.rows != 1:
                required_dimensions.add((op.name, 0))
            if op.columns != 1:
                required_dimensions.add((op.name, 1))

    # print(required_dimensions)

    i = 1
    for dims in dependent_dimensions(equations):
        if not dims.isdisjoint(required_dimensions):
            # print(dims)
            var = "n{}".format(i)
            i += 1
            for op_name, dim in dims:
                dimensions_dict[(op_name, dim)] = var

    # print(dimensions_dict)
    var_values = dict() # mapping of variable names to their values
    op_declarations = []
    const_names = dict()
    for expr in operands:
        props = [p.value for p in expr.properties if p.value in _SUPPORTED_PROPERTIES]
        prop_str = ', '.join(props)
        if isinstance(expr, ae.Matrix):
            rows = dimensions_dict[(expr.name, 0)]
            columns = dimensions_dict[(expr.name, 1)]
            var_values[rows] = expr.rows
            var_values[columns] = expr.columns

            # TODO for identity and zero matrix, we cannot use their original name.
            # Use something like I_n2
            if expr.has_property(Property.CONSTANT):
                if rows != columns:
                    name_suffix = "_{}_{}".format(rows, columns)
                else:
                    name_suffix = "_{}".format(rows)

            if expr.has_property(Property.IDENTITY):
                op_declarations.append(
                    'IdentityMatrix I{0}({1}, {2})'.format(name_suffix, rows, columns))
                const_names[expr] = "I{}".format(name_suffix)
            elif expr.has_property(Property.ZERO):
                op_declarations.append(
                    'ZeroMatrix Zero{0}({1}, {2})'.format(name_suffix, rows, columns))
                const_names[expr] = "Zero{}".format(name_suffix)
            else:
                op_declarations.append(
                    'Matrix {0}({1}, {2}) <{3}>'.format(expr.name, rows, columns, prop_str))
        elif isinstance(expr, ae.Vector):
            if expr.columns == 1:
                length = dimensions_dict[(expr.name, 0)]
                var_values[length] = expr.rows
                op_declarations.append('ColumnVector {0}({1}) <{2}>'.format(
                    expr.name, dimensions_dict[(expr.name, 0)], prop_str))
            else:
                length = dimensions_dict[(expr.name, 1)]
                var_values[length] = expr.columns
                op_declarations.append('RowVector {0}({1}) <{2}>'.format(
                    expr.name, dimensions_dict[(expr.name, 1)], prop_str))
        elif isinstance(expr, ae.Scalar):
            op_declarations.append('Scalar {} <{}>'.format(expr.name, prop_str))
        else:
            raise TypeError(
                'Unsupported operand type: {!r}'.format(type(expr)))

    var_declarations = []
    for variable, value in var_values.items():
        var_declarations.append('{} = {}'.format(variable, value))

    assignments = [export_expression(e, const_names) for e in equations]

    return '{}\n\n{}\n\n{}'.format(
        '\n'.join(sorted(var_declarations)),
        '\n'.join(sorted(op_declarations)),
        '\n'.join(assignments)
    )


def export_expression(expression, const_names):
    """Convert an expression to the input string format."""
    # TODO: Add support for conjugate etc. once they are in the input language
    if isinstance(expression, ae.Symbol):
        if expression in const_names:
            return const_names[expression]
        else:
            return expression.name
    if isinstance(expression, ae.ConstantScalar):
        return str(expression.value)
    if isinstance(expression, ae.Equal):
        return '{} = {}'.format(
            export_expression(expression.lhs, const_names),
            export_expression(expression.rhs, const_names)
        )
    for op_type, op_str in _BINARY.items():
        if isinstance(expression, op_type):
            return '({})'.format(op_str.join(export_expression(e, const_names) for e in expression.operands))
    for op_type, op_str in _UNARY.items():
        if isinstance(expression, op_type):
            return '{}({})'.format(op_str, export_expression(expression.operand, const_names))
    if isinstance(expression, ae.InverseTranspose):
        return 'inv(trans({}))'.format(export_expression(expression.operand, const_names))
    raise TypeError("Unknown expression type: {!r}".format(type(expression)))
