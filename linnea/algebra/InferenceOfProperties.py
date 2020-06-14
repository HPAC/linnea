import copy

from . import expression as ae
from .properties import Property, implications, negative_implications
from .. import temporaries

def isAuxiliary(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.AUXILIARY, isAuxiliary)
    return False

def isIdentity(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.IDENTITY, isIdentity)
    if isinstance(expr, ae.Plus):
        return False
    if isinstance(expr, ae.Times):
        return all(map(isIdentity, expr.operands))
    if isinstance(expr, (ae.Transpose, ae.Inverse, ae.InverseTranspose)):
        return isIdentity(expr.operand)
    return False

def isConstant(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.CONSTANT, isConstant)
    if isinstance(expr, ae.Operator):
        return all(map(isConstant, expr.operands))

def isFactor(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.FACTOR, isFactor)
    if isinstance(expr, (ae.Transpose, ae.Inverse, ae.InverseTranspose)):
        return isFactor(expr.operand)
    return False

def isUnitDiagonal(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.UNIT_DIAGONAL, isUnitDiagonal)
    if isinstance(expr, ae.Plus):
        return False
    if isinstance(expr, ae.Times): # TODO Should check triangular as well?
        return all(map(isUnitDiagonal, expr.operands)) and \
                (all(map(isLowerTriangular, expr.operands)) or
                 all(map(isUpperTriangular, expr.operands)))
    if isinstance(expr, (ae.Transpose, ae.Inverse, ae.InverseTranspose)):
        return isUnitDiagonal(expr.operand)
    return False

def isPositive(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.POSITIVE, isPositive)
    if isinstance(expr, (ae.Plus, ae.Times)):
        return all(map(isPositive, expr.operands))
    if isinstance(expr, (ae.Transpose, ae.Inverse, ae.InverseTranspose)):
        return isPositive(expr.operand)
    return False

def isSymmetric(expr):
    if isinstance(expr, ae.Symbol):
        if infer_property_symbol(expr, Property.SYMMETRIC, isSymmetric):
            return True
        # if expr is not symmetric (stored in false_properties), the two test below are still executed. Can this be avoided? Is it necessary? 
        elif isSquare(expr) and isDiagonalB(expr):
            expr.properties.add(Property.SYMMETRIC)
            return True
        else:
            return False
    else:
        # TODO As a shortcut, test if square, test bandwidth?
        return expr.transpose_of(expr)


def isSPSD(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.SPSD, isSPSD)
    if isinstance(expr, ae.Plus):
        return all(map(isSPSD, expr.operands))
    if isinstance(expr, ae.Times):
        return isSPSDTimes(expr)
    if isinstance(expr, (ae.Transpose, ae.Inverse, ae.InverseTranspose)):
        return isSPSD(expr.operand)
    return False


def isSPSDTimes(expr):
    scalars, non_scalars = expr.split_operands()
    if scalars and not isPositive(ae.Times(*scalars)):
        return False
    length = len(non_scalars)
    if length == 0:
        return False
    elif length == 1:
        return isSPSD(non_scalars[0])
    elif length % 2 == 0:
        left = ae.Times(*non_scalars[:length//2])
        right = ae.Times(*non_scalars[length//2:])
        return left.transpose_of(right)
    else:
        left = ae.Times(*non_scalars[:length//2])
        right = ae.Times(*non_scalars[length//2+1:])
        middle = non_scalars[length//2]
        return isSPSD(middle) and left.transpose_of(right)


def isSPD(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.SPD, isSPD)
    if isinstance(expr, ae.Plus):
        return all(map(isSPD, expr.operands))
    if isinstance(expr, ae.Times): # related to "iif they commute" ... ?
        return isSPDTimes(expr)
    if isinstance(expr, (ae.Transpose, ae.Inverse, ae.InverseTranspose)):
        return isSPD(expr.operand)
    return False


def isSPDTimes(expr):
    scalars, non_scalars = expr.split_operands()
    if scalars and not isPositive(ae.Times(*scalars)):
        return False
    length = len(non_scalars)
    if length == 0:
        return False
    elif length == 1:
        return isSPD(non_scalars[0])
    elif length % 2 == 0:
        left = ae.Times(*non_scalars[:length//2])
        right = ae.Times(*non_scalars[length//2:])
        return left.columns >= left.rows and isFullRank(left) and left.transpose_of(right)
    else:
        left = ae.Times(*non_scalars[:length//2])
        right = ae.Times(*non_scalars[length//2+1:])
        middle = non_scalars[length//2]
        return left.columns >= left.rows and isFullRank(left) and isSPD(middle) and left.transpose_of(right)


def isNonSingular(expr):
    if isinstance(expr, ae.Symbol):
        # if infer_property_symbol(expr, Property.SQUARE, isSquare) and infer_property_symbol(expr, Property.FULL_RANK, isFullRank):
        if isSquare(expr) and isFullRank(expr):
            expr.properties.add(Property.NON_SINGULAR)
            return True
        else:
            return infer_property_symbol(expr, Property.NON_SINGULAR, isNonSingular)
    if isinstance(expr, ae.Times):
        return all(map(isNonSingular, expr.operands))
    if isinstance(expr, ae.Plus): # ?
        return any(map(isNonSingular, expr.operands))
    if isinstance(expr, (ae.Transpose, ae.Inverse, ae.InverseTranspose)):
        return isNonSingular(expr.operand)
    return False

def isOrthogonal(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.ORTHOGONAL, isOrthogonal)
    if isinstance(expr, ae.Times):
        return all(map(isOrthogonal, expr.operands))
    if isinstance(expr, (ae.Transpose, ae.Inverse, ae.InverseTranspose)):
        return isOrthogonal(expr.operand)
    return False

def isOrthogonalColumns(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.ORTHOGONAL_COLUMNS, isOrthogonalColumns)
    if isinstance(expr, ae.Transpose):
        return isOrthogonalColumns(expr.operand)
    return False

def isOrthogonalRows(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.ORTHOGONAL_ROWS, isOrthogonalRows)
    if isinstance(expr, ae.Transpose):
        return isOrthogonalRows(expr.operand)
    return False

def isFullRank(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.FULL_RANK, isFullRank)
    if isinstance(expr, ae.Times):
        return all(map(isFullRank, expr.operands)) and is_full_rank_product(expr)
    if isinstance(expr, ae.Plus):
        return any(map(isFullRank, expr.operands))
    if isinstance(expr, (ae.Transpose, ae.Inverse, ae.InverseTranspose)):
        return isFullRank(expr.operand)
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

def isSquare(expr):
    return infer_property_test_function(expr, Property.SQUARE, isSquareTF)

def isSquareTF(expr):
    if isinstance(expr, ae.Scalar):
        return False
    else:
        size = expr.size
        if size[0] == size[1]:
            return True
        else:
            return False

def isColumnPanel(expr):
    return infer_property_test_function(expr, Property.COLUMN_PANEL, isColumnPanelTF)

def isColumnPanelTF(expr):
    if isinstance(expr, ae.Scalar):
        return False
    else:
        size = expr.size
        if size[0] > size[1]:
            return True
        else:
            return False

def isRowPanel(expr):
    return infer_property_test_function(expr, Property.ROW_PANEL, isRowPanelTF)

def isRowPanelTF(expr):
    if isinstance(expr, ae.Scalar):
        return False
    else:
        size = expr.size
        if size[0] < size[1]:
            return True
        else:
            return False

def isUnitary(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.UNITARY, isUnitary)
    return False
        
def isNormal(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.NORMAL, isNormal)
    return False
        
def isHermitian(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.HERMITIAN, isHermitian)
    return False
        
def isPermutation(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, Property.PERMUTATION, isPermutation)
    if isinstance(expr, ae.Times):
        return all(map(isPermutation, expr.operands))
    if isinstance(expr, (ae.Transpose, ae.Inverse, ae.InverseTranspose)):
        return isPermutation(expr.operand)
    return False

def admitsFactorization(expr):
    return infer_property_test_function(expr, Property.ADMITS_FACTORIZATION, admitsFactorizationTF)

def admitsFactorizationTF(expr):
    return not (isDiagonalB(expr) or isTriangularB(expr) or isOrthogonal(expr) or isOrthogonalColumns(expr) or isOrthogonalRows(expr) or isPermutation(expr))

def isScalar(expr):
    return infer_property_test_function(expr, Property.SCALAR, isScalarTF)

def isScalarTF(expr):
    size = expr.size
    if size[0] == 1 and size[1] == 1:
        return True
    else:
        return False

def isVector(expr):
    return infer_property_test_function(expr, Property.VECTOR, isVectorTF)

def isVectorTF(expr):
    rows, columns = expr.size
    if (rows == 1 and columns != 1) or (rows != 1 and columns == 1):
        return True
    else:
        return False

def isMatrix(expr):
    return infer_property_test_function(expr, Property.MATRIX, isMatrixTF)

def isMatrixTF(expr):
    size = expr.size
    if size[0] != 1 and size[1] != 1:
        return True
    else:
        return False

def isLowerTriangular(expr):
    return infer_property_test_function(expr, Property.LOWER_TRIANGULAR, isLowerTriangularB)

def isLowerTriangularB(expr):
    lb, ub = expr.bandwidth
    if ub <= 0:
        return True
    else:
        return False

def isUpperTriangular(expr):
    return infer_property_test_function(expr, Property.UPPER_TRIANGULAR, isUpperTriangularB)

def isUpperTriangularB(expr):
    lb, ub = expr.bandwidth
    if lb <= 0:
        return True
    else:
        return False

def isTriangular(expr):
    return infer_property_test_function(expr, Property.TRIANGULAR, isTriangularB)

def isTriangularB(expr):
    lb, ub = expr.bandwidth
    if lb <= 0 or ub <= 0:
        return True
    else:
        return False

def isDiagonal(expr):
    return infer_property_test_function(expr, Property.DIAGONAL, isDiagonalB)

def isDiagonalB(expr):
    lb, ub = expr.bandwidth
    if lb == 0 and ub == 0:
        return True
    else:
        return False

def isZero(expr):
    return infer_property_test_function(expr, Property.ZERO, isZeroB)

def isZeroB(expr):
    # Careful. Unless the bandwidth of the Zero matrix is set to something
    # appropriate (e.g. (0, -1)), this is not a property that purely depends
    # on bandwidth.
    lb, ub = expr.bandwidth
    if lb + ub + 1 <= 0:
        return True
    else:
        return False


def infer_property_test_function(expr, prop, test_func):
    if isinstance(expr, ae.Symbol):
        properties = expr.properties
        false_properties = expr.false_properties
        if prop in properties:
            return True
        if prop in false_properties:
            return False
        if test_func(expr):
            properties.add(prop)
            return True
        else:
            false_properties.add(prop)
            return False
    else:
        return test_func(expr)


def infer_property_symbol(expr, prop, test_func):
    """
    TODO: Use partial ordering directly? Is that necessary if the implications
    are used (they are based on the ordering)?
    """
    properties = expr.properties
    false_properties = expr.false_properties
    if prop in properties:
        return True
    elif prop in false_properties:
        return False
    else:
        try:
            equivalent_expr = temporaries._equivalent_expressions[expr.name]
        except KeyError:
            return False
        else:
            # print(expr, equivalent_expr, prop, infer_property(equivalent_expr, prop))
            has_property = test_func(equivalent_expr)
            if has_property:
                properties.add(prop)
                properties.update(implications.get(prop, tuple()))
                false_properties.update(negative_implications.get(prop, tuple()))
            else:
                false_properties.add(prop)
            return has_property


property_to_function = {
    Property.AUXILIARY: isAuxiliary,
    Property.ZERO: isZero,
    Property.IDENTITY: isIdentity,
    Property.DIAGONAL: isDiagonal,
    Property.TRIANGULAR: isTriangular,
    Property.LOWER_TRIANGULAR: isLowerTriangular,
    Property.UPPER_TRIANGULAR: isUpperTriangular,
    Property.UNIT_DIAGONAL: isUnitDiagonal,
    Property.POSITIVE: isPositive,
    Property.SYMMETRIC: isSymmetric,
    Property.SPSD: isSPSD,
    Property.SPD: isSPD,
    Property.NON_SINGULAR: isNonSingular,
    Property.ORTHOGONAL: isOrthogonal,
    Property.FULL_RANK: isFullRank,
    Property.SQUARE: isSquare,
    Property.SCALAR: isScalar,
    Property.VECTOR: isVector,
    Property.MATRIX: isMatrix,
    Property.COLUMN_PANEL: isColumnPanel,
    Property.ROW_PANEL: isRowPanel,
    Property.UNITARY: isUnitary,
    Property.NORMAL: isNormal,
    Property.HERMITIAN: isHermitian,
    Property.ORTHOGONAL_COLUMNS: isOrthogonalColumns,
    Property.ORTHOGONAL_ROWS: isOrthogonalRows,
    Property.PERMUTATION: isPermutation,
    Property.ADMITS_FACTORIZATION: admitsFactorization,
    Property.FACTOR: isFactor,
    Property.CONSTANT: isConstant,
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