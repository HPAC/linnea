import copy

from . import expression as ae
from .properties import Property, implications, negative_implications, \
                        binary_implications_backwards
from .. import temporaries

def is_auxiliary(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.AUXILIARY, is_auxiliary)
    return False

def is_identity(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.IDENTITY, is_identity)
    if isinstance(expr, ae.Plus):
        return False
    if isinstance(expr, ae.Times):
        return all(map(is_identity, expr.operands))
    if isinstance(expr, (ae.Transpose, ae.Inverse, ae.InverseTranspose)):
        return is_identity(expr.operand)
    return False

def is_constant(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.CONSTANT, is_constant)
    if isinstance(expr, ae.Operator):
        return all(map(is_constant, expr.operands))

def is_factor(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.FACTOR, is_factor)
    if isinstance(expr, (ae.Transpose, ae.Inverse, ae.InverseTranspose)):
        return is_factor(expr.operand)
    return False

def is_unit_diagonal(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.UNIT_DIAGONAL, is_unit_diagonal)
    if isinstance(expr, ae.Plus):
        return False
    if isinstance(expr, ae.Times): # TODO Should check triangular as well?
        return all(map(is_unit_diagonal, expr.operands)) and \
                (all(map(is_lower_triangular, expr.operands)) or
                 all(map(is_upper_triangular, expr.operands)))
    if isinstance(expr, (ae.Transpose, ae.Inverse, ae.InverseTranspose)):
        return is_unit_diagonal(expr.operand)
    return False

def is_positive(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.POSITIVE, is_positive)
    if isinstance(expr, (ae.Plus, ae.Times)):
        return all(map(is_positive, expr.operands))
    if isinstance(expr, (ae.Transpose, ae.Inverse, ae.InverseTranspose)):
        return is_positive(expr.operand)
    return False

def is_symmetric(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.SYMMETRIC, is_symmetric)
    else:
        # TODO As a shortcut, test if square, test bandwidth?
        return expr.transpose_of(expr)


def is_SPSD(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.SPSD, is_SPSD)
    if isinstance(expr, ae.Plus):
        return all(map(is_SPSD, expr.operands))
    if isinstance(expr, ae.Times):
        return is_SPSD_product(expr)
    if isinstance(expr, (ae.Transpose, ae.Inverse, ae.InverseTranspose)):
        return is_SPSD(expr.operand)
    return False


def is_SPSD_product(expr):
    scalars, non_scalars = expr.split_operands()
    if scalars and not is_positive(ae.Times(*scalars)):
        return False
    length = len(non_scalars)
    if length == 0:
        return False
    elif length == 1:
        return is_SPSD(non_scalars[0])
    elif length % 2 == 0:
        left = ae.Times(*non_scalars[:length//2])
        right = ae.Times(*non_scalars[length//2:])
        return left.transpose_of(right)
    else:
        left = ae.Times(*non_scalars[:length//2])
        right = ae.Times(*non_scalars[length//2+1:])
        middle = non_scalars[length//2]
        return is_SPSD(middle) and left.transpose_of(right)


def is_SPD(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.SPD, is_SPD)
    if isinstance(expr, ae.Plus):
        return all(map(is_SPD, expr.operands))
    if isinstance(expr, ae.Times): # related to "iif they commute" ... ?
        return is_SPD_product(expr)
    if isinstance(expr, (ae.Transpose, ae.Inverse, ae.InverseTranspose)):
        return is_SPD(expr.operand)
    return False


def is_SPD_product(expr):
    scalars, non_scalars = expr.split_operands()
    if scalars and not is_positive(ae.Times(*scalars)):
        return False
    length = len(non_scalars)
    if length == 0:
        return False
    elif length == 1:
        return is_SPD(non_scalars[0])
    elif length % 2 == 0:
        left = ae.Times(*non_scalars[:length//2])
        right = ae.Times(*non_scalars[length//2:])
        return left.columns >= left.rows and is_full_rank(left) and left.transpose_of(right)
    else:
        left = ae.Times(*non_scalars[:length//2])
        right = ae.Times(*non_scalars[length//2+1:])
        middle = non_scalars[length//2]
        return left.columns >= left.rows and is_full_rank(left) and is_SPD(middle) and left.transpose_of(right)


def is_nonsingular(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.NON_SINGULAR, is_nonsingular)
    if isinstance(expr, ae.Times):
        return all(map(is_nonsingular, expr.operands))
    if isinstance(expr, ae.Plus): # ?
        return any(map(is_nonsingular, expr.operands))
    if isinstance(expr, (ae.Transpose, ae.Inverse, ae.InverseTranspose)):
        return is_nonsingular(expr.operand)
    return False

def is_orthogonal(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.ORTHOGONAL, is_orthogonal)
    if isinstance(expr, ae.Times):
        return all(map(is_orthogonal, expr.operands))
    if isinstance(expr, (ae.Transpose, ae.Inverse, ae.InverseTranspose)):
        return is_orthogonal(expr.operand)
    return False

def is_orthogonal_columns(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.ORTHOGONAL_COLUMNS, is_orthogonal_columns)
    if isinstance(expr, ae.Transpose):
        return is_orthogonal_columns(expr.operand)
    return False

def is_orthogonal_rows(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.ORTHOGONAL_ROWS, is_orthogonal_rows)
    if isinstance(expr, ae.Transpose):
        return is_orthogonal_rows(expr.operand)
    return False

def is_full_rank(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.FULL_RANK, is_full_rank)
    if isinstance(expr, ae.Times):
        return all(map(is_full_rank, expr.operands)) and is_full_rank_product(expr)
    if isinstance(expr, ae.Plus):
        return any(map(is_full_rank, expr.operands))
    if isinstance(expr, (ae.Transpose, ae.Inverse, ae.InverseTranspose)):
        return is_full_rank(expr.operand)
    return False

def is_full_rank_product(expr):
    """Tests if product is full rank based on operand sizes.

    Args:
        expr (Expression): This expression is assumed to be a product.

    Returns:
        True if product is full rank, False if not.

    TODO:
        - This function does not work correctly with inner products.
        - Consider bandwidth?
    """
    max_rank = min(expr.size)
    _, non_scalars = expr.split_operands()
    min_interior_size = min(operand.columns for operand in non_scalars[:-1])
    return min_interior_size >= max_rank

def is_square(expr):
    return infer_property_test_function(expr, Property.SQUARE, is_square_TF)

def is_square_TF(expr):
    if isinstance(expr, ae.Scalar):
        return False
    else:
        size = expr.size
        if size[0] == size[1]:
            return True
        else:
            return False

def is_columnpanel(expr):
    return infer_property_test_function(expr, Property.COLUMN_PANEL, is_columnpanel_TF)

def is_columnpanel_TF(expr):
    if isinstance(expr, ae.Scalar):
        return False
    else:
        size = expr.size
        if size[0] > size[1]:
            return True
        else:
            return False

def is_rowpanel(expr):
    return infer_property_test_function(expr, Property.ROW_PANEL, is_rowpanel_TF)

def is_rowpanel_TF(expr):
    if isinstance(expr, ae.Scalar):
        return False
    else:
        size = expr.size
        if size[0] < size[1]:
            return True
        else:
            return False

def is_unitary(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.UNITARY, is_unitary)
    return False
        
def is_normal(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.NORMAL, is_normal)
    return False
        
def is_hermitian(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.HERMITIAN, is_hermitian)
    return False
        
def is_permutation(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.PERMUTATION, is_permutation)
    if isinstance(expr, ae.Times):
        return all(map(is_permutation, expr.operands))
    if isinstance(expr, (ae.Transpose, ae.Inverse, ae.InverseTranspose)):
        return is_permutation(expr.operand)
    return False

def admits_factorization(expr):
    return infer_property_test_function(expr, Property.ADMITS_FACTORIZATION, admits_factorization_TF)

def admits_factorization_TF(expr):
    """Tests if an factorization can be applied to an expression.

    Factorizations can only be applied to matrices that are not
    - upper or lower triangular
    - diagonal
    - permutation matrices
    - orthogonal
    - have orthogonal rows or columns

    Since diagonal is a special case of triangular, and both orthogonal and
    permutation matrices have orthogonal rows and columns, it is sufficient to
    check that the expression is not triangular or has orthogonal rows or
    columns.
    """
    if not is_matrix(expr):
        return False
    return not (is_triangular_B(expr) or is_orthogonal_columns(expr) or is_orthogonal_rows(expr))

def is_scalar(expr):
    return infer_property_test_function(expr, Property.SCALAR, is_scalar_TF)

def is_scalar_TF(expr):
    size = expr.size
    if size[0] == 1 and size[1] == 1:
        return True
    else:
        return False

def is_vector(expr):
    return infer_property_test_function(expr, Property.VECTOR, is_vector_TF)

def is_vector_TF(expr):
    rows, columns = expr.size
    if (rows == 1 and columns != 1) or (rows != 1 and columns == 1):
        return True
    else:
        return False

def is_matrix(expr):
    return infer_property_test_function(expr, Property.MATRIX, is_matrix_TF)

def is_matrix_TF(expr):
    size = expr.size
    if size[0] != 1 and size[1] != 1:
        return True
    else:
        return False

def is_lower_triangular(expr):
    return infer_property_test_function(expr, Property.LOWER_TRIANGULAR, is_lower_triangular_B)

def is_lower_triangular_B(expr):
    lb, ub = expr.bandwidth
    if ub <= 0:
        return True
    else:
        return False

def is_upper_triangular(expr):
    return infer_property_test_function(expr, Property.UPPER_TRIANGULAR, is_upper_triangular_B)

def is_upper_triangular_B(expr):
    lb, ub = expr.bandwidth
    if lb <= 0:
        return True
    else:
        return False

def is_triangular(expr):
    return infer_property_test_function(expr, Property.TRIANGULAR, is_triangular_B)

def is_triangular_B(expr):
    lb, ub = expr.bandwidth
    if lb <= 0 or ub <= 0:
        return True
    else:
        return False

def is_diagonal(expr):
    return infer_property_test_function(expr, Property.DIAGONAL, is_diagonal_B)

def is_diagonal_B(expr):
    lb, ub = expr.bandwidth
    if lb == 0 and ub == 0:
        return True
    else:
        return False

def is_zero(expr):
    return infer_property_test_function(expr, Property.ZERO, is_zero_B)

def is_zero_B(expr):
    # Careful. Unless the bandwidth of the Zero matrix is set to something
    # appropriate (e.g. (0, -1)), this is not a property that purely depends
    # on bandwidth.
    lb, ub = expr.bandwidth
    if lb + ub + 1 <= 0:
        return True
    else:
        return False


def add_property(expr, property, has_property):
    if has_property:
        expr.properties.add(property)
        expr.properties.update(implications.get(property, set()))
        expr.false_properties.update(negative_implications.get(property, set()))
    else:
        expr.false_properties.add(property)


def infer_property_test_function(expr, prop, test_func):
    if isinstance(expr, ae.Symbol):
        if prop in expr.properties:
            return True
        if prop in expr.false_properties:
            return False
        has_property = test_func(expr)
        add_property(expr, prop, has_property)
        return has_property
    else:
        return test_func(expr)


def infer_property_symbol(expr, prop, test_func):
    if prop in expr.properties:
        return True
    elif prop in expr.false_properties:
        return False

    try:
        equivalent_expr = temporaries._equivalent_expressions[expr.name]
    except KeyError:
        pass
    else:
        # print(expr, equivalent_expr, prop, infer_property(equivalent_expr, prop))
        has_property = test_func(equivalent_expr)
        add_property(expr, prop, has_property)
        if has_property:
            return True

    try:
        property_sets = binary_implications_backwards[prop]
    except KeyError:
        pass
    else:
        for required_properties in property_sets:
            # print(prop, ss)
            has_property = all(property_to_function[prop](expr) for prop in required_properties)
            # print(r)
            add_property(expr, prop, has_property)
            # If is necessary here because of the loop. There could be more than
            # one rule.
            if has_property:
                return True
    return False
            


property_to_function = {
    Property.AUXILIARY: is_auxiliary,
    Property.ZERO: is_zero,
    Property.IDENTITY: is_identity,
    Property.DIAGONAL: is_diagonal,
    Property.TRIANGULAR: is_triangular,
    Property.LOWER_TRIANGULAR: is_lower_triangular,
    Property.UPPER_TRIANGULAR: is_upper_triangular,
    Property.UNIT_DIAGONAL: is_unit_diagonal,
    Property.POSITIVE: is_positive,
    Property.SYMMETRIC: is_symmetric,
    Property.SPSD: is_SPSD,
    Property.SPD: is_SPD,
    Property.NON_SINGULAR: is_nonsingular,
    Property.ORTHOGONAL: is_orthogonal,
    Property.FULL_RANK: is_full_rank,
    Property.SQUARE: is_square,
    Property.SCALAR: is_scalar,
    Property.VECTOR: is_vector,
    Property.MATRIX: is_matrix,
    Property.COLUMN_PANEL: is_columnpanel,
    Property.ROW_PANEL: is_rowpanel,
    Property.UNITARY: is_unitary,
    Property.NORMAL: is_normal,
    Property.HERMITIAN: is_hermitian,
    Property.ORTHOGONAL_COLUMNS: is_orthogonal_columns,
    Property.ORTHOGONAL_ROWS: is_orthogonal_rows,
    Property.PERMUTATION: is_permutation,
    Property.ADMITS_FACTORIZATION: admits_factorization,
    Property.FACTOR: is_factor,
    Property.CONSTANT: is_constant,
}

def infer_property(expr, prop):
    if not expr.is_constant: # does the expression contain variables?
        return False
    try:
        func = property_to_function[prop]
    except KeyError:
        raise NotImplementedError("No function to infer %s" % prop)
    return func(expr)

if __name__ == "__main__":
    pass