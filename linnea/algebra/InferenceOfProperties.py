import copy

from . import expression as ae
from .properties import Property as properties
from .properties import implications, negative_implications
from .. import temporaries

def isInput(expr):
    # isinstance?
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, properties.INPUT, isInput)
    if isinstance(expr, ae.Plus):
        return all(isInput(term) for term in expr.operands)
    if isinstance(expr, ae.Times):
        return all(isInput(factor) for factor in expr.operands)
    if isinstance(expr, ae.Transpose):
        return isInput(expr.operand)
    if isinstance(expr, ae.Inverse):
        return isInput(expr.operand)
    if isinstance(expr, ae.InverseTranspose):
        return isInput(expr.operand)
    # [TODO] Double check!!!
    return False


def isOutput(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, properties.INPUT, isOutput)
    if isinstance(expr, ae.Plus):
        return any(isOutput(term) for term in expr.operands)
    if isinstance(expr, ae.Times):
        return any(isOutput(factor) for factor in expr.operands)
    if isinstance(expr, ae.Transpose):
        return isOutput(expr.operand)
    if isinstance(expr, ae.Inverse):
        return isOutput(expr.operand)
    if isinstance(expr, ae.InverseTranspose):
        return isOutput(expr.operand)
    return False

def isAuxiliary(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, properties.AUXILIARY, isAuxiliary)
    return False

def isIdentity(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, properties.IDENTITY, isIdentity)
    if isinstance(expr, ae.Plus):
        return False
    if isinstance(expr, ae.Times):
        return all(isIdentity(factor) for factor in expr.operands)
    if isinstance(expr, ae.Transpose):
        return isIdentity(expr.operand)
    if isinstance(expr, ae.Inverse):
        return isIdentity(expr.operand)
    if isinstance(expr, ae.InverseTranspose):
        return isIdentity(expr.operand)
    return False

def isConstant(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, properties.CONSTANT, isConstant)
    if isinstance(expr, ae.Operator):
        return all(isConstant(operand) for operand in expr.operands)

def isFactor(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, properties.FACTOR, isFactor)
    if isinstance(expr, ae.Transpose):
        return isFactor(expr.operand)
    if isinstance(expr, ae.Inverse):
        return isFactor(expr.operand)
    if isinstance(expr, ae.InverseTranspose):
        return isFactor(expr.operand)
    return False

def isUnitDiagonal(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, properties.UNIT_DIAGONAL, isUnitDiagonal)
    if isinstance(expr, ae.Plus):
        return False
    if isinstance(expr, ae.Times): # TODO Should check triangular as well?
        return all(isUnitDiagonal(factor) for factor in expr.operands) and \
                (all(isLowerTriangular(factor) for factor in expr.operands) or
                 all(isUpperTriangular(factor) for factor in expr.operands))
    if isinstance(expr, ae.Transpose):
        return isUnitDiagonal(expr.operand)
    if isinstance(expr, ae.Inverse): # TODO triangular?
        return isUnitDiagonal(expr.operand)
    if isinstance(expr, ae.InverseTranspose): # TODO triangular?
        return isUnitDiagonal(expr.operand)
    return False

def isPositive(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, properties.POSITIVE, isPositive)
    if isinstance(expr, ae.Plus):
        return all(isPositive(term) for term in expr.operands)
    if isinstance(expr, ae.Times):
        return all(isPositive(term) for term in expr.operands)
    if isinstance(expr, ae.Transpose):
        return isPositive(expr.operand)
    if isinstance(expr, ae.Inverse):
        return isPositive(expr.operand)
    if isinstance(expr, ae.InverseTranspose):
        return isPositive(expr.operand)
    return False

def isSymmetric(expr):
    if isinstance(expr, ae.Symbol):
        if infer_property_symbol(expr, properties.SYMMETRIC, isSymmetric):
            return True
        # if expr is not symmetric (stored in false_properties), the two test below are still executed. Can this be avoided? Is it necessary? 
        elif isSquare(expr) and isDiagonalB(expr):
            expr.properties.add(properties.SYMMETRIC)
            return True
        else:
            return False
    else:
        # TODO As a shortcut, test if square, test bandwidth?
        return expr.transpose_of(expr)


def isSPSD(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, properties.SPSD, isSPSD)
    if isinstance(expr, ae.Plus):
        return all(isSPSD(term) for term in expr.operands)
    if isinstance(expr, ae.Times):
        return isSPSDTimes(expr)
    if isinstance(expr, ae.Transpose):
        return isSPSD(expr.operand)
    if isinstance(expr, ae.Inverse):
        return isSPSD(expr.operand)
    if isinstance(expr, ae.InverseTranspose):
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
        return infer_property_symbol(expr, properties.SPD, isSPD)
    if isinstance(expr, ae.Plus):
        return all(isSPD(term) for term in expr.operands)
    if isinstance(expr, ae.Times): # related to "iif they commute" ... ?
        return isSPDTimes(expr)
    if isinstance(expr, ae.Transpose):
        return isSPD(expr.operand)
    if isinstance(expr, ae.Inverse):
        return isSPD(expr.operand)
    if isinstance(expr, ae.InverseTranspose):
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
        # if infer_property_symbol(expr, properties.SQUARE, isSquare) and infer_property_symbol(expr, properties.FULL_RANK, isFullRank):
        if isSquare(expr) and isFullRank(expr):
            expr.properties.add(properties.NON_SINGULAR)
            return True
        else:
            return infer_property_symbol(expr, properties.NON_SINGULAR, isNonSingular)
    if isinstance(expr, ae.Times):
        return all(isNonSingular(factor) for factor in expr.operands)
    if isinstance(expr, ae.Plus): # ?
        return all(isNonSingular(operand) for operand in expr.operands)
    if isinstance(expr, ae.Transpose):
        return isNonSingular(expr.operand)
    if isinstance(expr, ae.Inverse):
        return isNonSingular(expr.operand)
    if isinstance(expr, ae.InverseTranspose):
        return isNonSingular(expr.operand)
    return False

def isOrthogonal(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, properties.ORTHOGONAL, isOrthogonal)
    if isinstance(expr, ae.Times):
        return all(isOrthogonal(factor) for factor in expr.operands)
    if isinstance(expr, ae.Transpose):
        return isOrthogonal(expr.operand)
    if isinstance(expr, ae.Inverse):
        return isOrthogonal(expr.operand)
    if isinstance(expr, ae.InverseTranspose):
        return isOrthogonal(expr.operand)
    return False

def isOrthogonalColumns(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, properties.ORTHOGONAL_COLUMNS, isOrthogonalColumns)
    if isinstance(expr, ae.Transpose):
        return isOrthogonalColumns(expr.operand)
    return False

def isOrthogonalRows(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, properties.ORTHOGONAL_ROWS, isOrthogonalRows)
    if isinstance(expr, ae.Transpose):
        return isOrthogonalRows(expr.operand)
    return False

def isFullRank(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, properties.FULL_RANK, isFullRank)
    if isinstance(expr, ae.Times):
        return all(isFullRank(factor) for factor in expr.operands) and is_full_rank_product(expr)
    if isinstance(expr, ae.Plus):
        return all(isFullRank(operand) for operand in expr.operands)
    if isinstance(expr, ae.Transpose):
        return isFullRank(expr.operand)
    if isinstance(expr, ae.Inverse):
        return isFullRank(expr.operand)
    if isinstance(expr, ae.InverseTranspose):
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
    min_interior_size = min(operand.size[1] for operand in non_scalars[:-1])
    return min_interior_size >= max_rank

def isSquare(expr):
    return infer_property_test_function(expr, properties.SQUARE, isSquareTF)

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
    return infer_property_test_function(expr, properties.COLUMN_PANEL, isColumnPanelTF)

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
    return infer_property_test_function(expr, properties.ROW_PANEL, isRowPanelTF)

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
        return infer_property_symbol(expr, properties.UNITARY, isUnitary)
    return False
        
def isNormal(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, properties.NORMAL, isNormal)
    return False
        
def isHermitian(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, properties.HERMITIAN, isHermitian)
    return False
        
def isPermutation(expr):
    if isinstance(expr, ae.Symbol):
        return infer_property_symbol(expr, properties.PERMUTATION, isPermutation)
    if isinstance(expr, ae.Times):
        return all(isPermutation(factor) for factor in expr.operands)
    if isinstance(expr, ae.Transpose):
        return isPermutation(expr.operand)
    if isinstance(expr, ae.Inverse):
        return isPermutation(expr.operand)
    if isinstance(expr, ae.InverseTranspose):
        return isPermutation(expr.operand)
    return False

def admitsFactorization(expr):
    return infer_property_test_function(expr, properties.ADMITS_FACTORIZATION, admitsFactorizationTF)

def admitsFactorizationTF(expr):
    return not (isDiagonalB(expr) or isTriangularB(expr) or isOrthogonal(expr) or isOrthogonalColumns(expr) or isOrthogonalRows(expr) or isPermutation(expr))

def isScalar(expr):
    return infer_property_test_function(expr, properties.SCALAR, isScalarTF)

def isScalarTF(expr):
    size = expr.size
    if size[0] == 1 and size[1] == 1:
        return True
    else:
        return False

def isVector(expr):
    return infer_property_test_function(expr, properties.VECTOR, isVectorTF)

def isVectorTF(expr):
    rows, columns = expr.size
    if (rows == 1 and columns != 1) or (rows != 1 and columns == 1):
        return True
    else:
        return False

def isMatrix(expr):
    return infer_property_test_function(expr, properties.MATRIX, isMatrixTF)

def isMatrixTF(expr):
    size = expr.size
    if size[0] != 1 and size[1] != 1:
        return True
    else:
        return False

def isLowerTriangular(expr):
    return infer_property_test_function(expr, properties.LOWER_TRIANGULAR, isLowerTriangularB)

def isLowerTriangularB(expr):
    lb, ub = expr.bandwidth
    if ub <= 0:
        return True
    else:
        return False

def isUpperTriangular(expr):
    return infer_property_test_function(expr, properties.UPPER_TRIANGULAR, isUpperTriangularB)

def isUpperTriangularB(expr):
    lb, ub = expr.bandwidth
    if lb <= 0:
        return True
    else:
        return False

def isTriangular(expr):
    return infer_property_test_function(expr, properties.TRIANGULAR, isTriangularB)

def isTriangularB(expr):
    lb, ub = expr.bandwidth
    if lb <= 0 or ub <= 0:
        return True
    else:
        return False

def isDiagonal(expr):
    return infer_property_test_function(expr, properties.DIAGONAL, isDiagonalB)

def isDiagonalB(expr):
    lb, ub = expr.bandwidth
    if lb == 0 and ub == 0:
        return True
    else:
        return False

def isZero(expr):
    return infer_property_test_function(expr, properties.ZERO, isZeroB)

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
        _properties = expr.properties
        _false_properties = expr.false_properties
        if prop in _properties:
            return True
        if prop in _false_properties:
            return False
        if test_func(expr):
            _properties.add(prop)
            return True
        else:
            _false_properties.add(prop)
            return False
    else:
        return test_func(expr)


def infer_property_symbol(expr, prop, test_func):
    """
    TODO: Use partial ordering directly? Is that necessary if the implications
    are used (they are based on the ordering)?
    """
    _properties = expr.properties
    _false_properties = expr.false_properties
    if prop in _properties:
        return True
    elif prop in _false_properties:
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
                _properties.add(prop)
                _properties.update(implications.get(prop, tuple()))
                _false_properties.update(negative_implications.get(prop, tuple()))
            else:
                _false_properties.add(prop)
            return has_property


property_to_function = {
    properties.INPUT: isInput,
    properties.OUTPUT: isOutput,
    properties.AUXILIARY: isAuxiliary,
    properties.ZERO: isZero,
    properties.IDENTITY: isIdentity,
    properties.DIAGONAL: isDiagonal,
    properties.TRIANGULAR: isTriangular,
    properties.LOWER_TRIANGULAR: isLowerTriangular,
    properties.UPPER_TRIANGULAR: isUpperTriangular,
    properties.UNIT_DIAGONAL: isUnitDiagonal,
    properties.POSITIVE: isPositive,
    properties.SYMMETRIC: isSymmetric,
    properties.SPSD: isSPSD,
    properties.SPD: isSPD,
    properties.NON_SINGULAR: isNonSingular,
    properties.ORTHOGONAL: isOrthogonal,
    properties.FULL_RANK: isFullRank,
    properties.SQUARE: isSquare,
    properties.SCALAR: isScalar,
    properties.VECTOR: isVector,
    properties.MATRIX: isMatrix,
    properties.COLUMN_PANEL: isColumnPanel,
    properties.ROW_PANEL: isRowPanel,
    properties.UNITARY: isUnitary,
    properties.NORMAL: isNormal,
    properties.HERMITIAN: isHermitian,
    properties.ORTHOGONAL_COLUMNS: isOrthogonalColumns,
    properties.ORTHOGONAL_ROWS: isOrthogonalRows,
    properties.PERMUTATION: isPermutation,
    properties.ADMITS_FACTORIZATION: admitsFactorization,
    properties.FACTOR: isFactor,
    properties.CONSTANT: isConstant,
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