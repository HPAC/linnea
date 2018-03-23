import copy

from . import expression as ae
from .properties import Property as properties
from .properties import implications, negative_implications
from .. import temporaries

from . import property_DNs

def isInput(node):
    # isinstance?
    if isinstance(node, ae.Symbol):
        return infer_property_symbol(node, properties.INPUT, isInput)
    if isinstance(node, ae.Plus):
        return all(isInput(term) for term in node.operands)
    if isinstance(node, ae.Times):
        return all(isInput(factor) for factor in node.operands)
    if isinstance(node, ae.Transpose):
        return isInput(node.operand)
    if isinstance(node, ae.Inverse):
        return isInput(node.operand)
    # [TODO] Double check!!!
    return False


def isOutput(node):
    if isinstance(node, ae.Symbol):
        return infer_property_symbol(node, properties.INPUT, isOutput)
    if isinstance(node, ae.Plus):
        return any(isOutput(term) for term in node.operands)
    if isinstance(node, ae.Times):
        return any(isOutput(factor) for factor in node.operands)
    if isinstance(node, ae.Transpose):
        return isOutput(node.operand)
    if isinstance(node, ae.Inverse):
        return isOutput(node.operand)
    return False

def isAuxiliary(node):
    if isinstance(node, ae.Symbol):
        return infer_property_symbol(node, properties.AUXILIARY, isAuxiliary)
    return False

def isIdentity(node):
    if isinstance(node, ae.Symbol):
        return infer_property_symbol(node, properties.IDENTITY, isIdentity)
    if isinstance(node, ae.Plus):
        return False
    if isinstance(node, ae.Times):
        return all(isIdentity(factor) for factor in node.operands)
    if isinstance(node, ae.Transpose):
        return isIdentity(node.operand)
    if isinstance(node, ae.Inverse):
        return isIdentity(node.operand)
    return False

def isConstant(node):
    if isinstance(node, ae.Symbol):
        return infer_property_symbol(node, properties.CONSTANT, isConstant)
    if isinstance(node, ae.Operator):
        return all(isConstant(operand) for operand in node.operands)

def isFactor(node):
    if isinstance(node, ae.Symbol):
        return infer_property_symbol(node, properties.FACTOR, isFactor)
    if isinstance(node, ae.Transpose):
        return isFactor(node.operand)
    if isinstance(node, ae.Inverse):
        return isFactor(node.operand)
    if isinstance(node, ae.InverseTranspose):
        return isFactor(node.operand)
    return False

def isUnitDiagonal(node):
    if isinstance(node, ae.Symbol):
        return infer_property_symbol(node, properties.UNIT_DIAGONAL, isUnitDiagonal)
    if isinstance(node, ae.Plus):
        return False
    if isinstance(node, ae.Times): # TODO Should check triangular as well?
        return all(isUnitDiagonal(factor) for factor in node.operands) and \
                (all(isLowerTriangular(factor) for factor in node.operands) or
                 all(isUpperTriangular(factor) for factor in node.operands))
    if isinstance(node, ae.Transpose):
        return isUnitDiagonal(node.operand)
    if isinstance(node, ae.Inverse): # TODO triangular?
        return isUnitDiagonal(node.operand)
    return False

def isSymmetric(node):
    if isinstance(node, ae.Symbol):
        if infer_property_symbol(node, properties.SYMMETRIC, isSymmetric):
            return True
        # if node is not symmetric (stored in false_properties), the two test below are still executed. Can this be avoided? Is it necessary? 
        elif isSquare(node) and isDiagonalB(node):
            node.properties.add(properties.SYMMETRIC)
            return True
        else:
            return False
    else:
        # TODO As a shortcut, test if square, test bandwidth?
        return node.transpose_of(node)


def isSPD(node):
    """
    TODO missing
    S + alpha I: If S is SPD and alpha is positive, the expression is SPD.
    What is missing: alpha I is currently not SPD. I should be SPD, for alpha,
    we need the property positive. A product is positive if all operands
    are positive, I is positive.
    Alterantively: In a sum, there is at least one SPD matrix and all the others
    are positive diagonal matrices. But positive and diagonal implies SPD.
    """
    # print(node)
    if isinstance(node, ae.Symbol):
        return infer_property_symbol(node, properties.SPD, isSPD)
    if isinstance(node, ae.Plus):
        return all(isSPD(term) for term in node.operands)
    if isinstance(node, ae.Times): # related to "iif they commute" ... ?
        return property_DNs.SPD_DN.is_match(node)
    if isinstance(node, ae.Transpose):
        return isSPD(node.operand)
    if isinstance(node, ae.Inverse):
        return isSPD(node.operand)
    return False

def isNonSingular(node):
    if isinstance(node, ae.Symbol):
        # if infer_property_symbol(node, properties.SQUARE, isSquare) and infer_property_symbol(node, properties.FULL_RANK, isFullRank):
        if isSquare(node) and isFullRank(node):
            node.properties.add(properties.NON_SINGULAR)
            return True
        else:
            return infer_property_symbol(node, properties.NON_SINGULAR, isNonSingular)
    if isinstance(node, ae.Plus): # ?
        return False
    if isinstance(node, ae.Times):
        return all(isNonSingular(factor) for factor in node.operands)
    if isinstance(node, ae.Transpose):
        return isNonSingular(node.operand)
    if isinstance(node, ae.Inverse):
        return isNonSingular(node.operand)
    return False

def isOrthogonal(node):
    if isinstance(node, ae.Symbol):
        return infer_property_symbol(node, properties.ORTHOGONAL, isOrthogonal)
    if isinstance(node, ae.Times):
        return all(isOrthogonal(factor) for factor in node.operands)
    if isinstance(node, ae.Transpose):
        return isOrthogonal(node.operand)
    if isinstance(node, ae.Inverse):
        return isOrthogonal(node.operand)
    return False

def isOrthogonalColumns(node):
    if isinstance(node, ae.Symbol):
        return infer_property_symbol(node, properties.ORTHOGONAL_COLUMNS, isOrthogonalColumns)
    if isinstance(node, ae.Transpose):
        return isOrthogonalColumns(node.operand)
    return False

def isOrthogonalRows(node):
    if isinstance(node, ae.Symbol):
        return infer_property_symbol(node, properties.ORTHOGONAL_ROWS, isOrthogonalRows)
    if isinstance(node, ae.Transpose):
        return isOrthogonalRows(node.operand)
    return False

def isFullRank(node):
    if isinstance(node, ae.Symbol):
        return infer_property_symbol(node, properties.FULL_RANK, isFullRank)
    if isinstance(node, ae.Times):
        return all(isFullRank(factor) for factor in node.operands) and is_full_rank_product(node)
    if isinstance(node, ae.Transpose):
        return isFullRank(node.operand)
    if isinstance(node, ae.Inverse):
        return isFullRank(node.operand)
    if isinstance(node, ae.InverseTranspose):
        return isFullRank(node.operand)
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

def isSquareTF(node):
    if isinstance(node, ae.Scalar):
        return False
    else:
        size = node.size
        if size[0] == size[1]:
            return True
        else:
            return False

def isColumnPanel(expr):
    return infer_property_test_function(expr, properties.COLUMN_PANEL, isColumnPanelTF)

def isColumnPanelTF(node):
    if isinstance(node, ae.Scalar):
        return False
    else:
        size = node.size
        if size[0] > size[1]:
            return True
        else:
            return False

def isRowPanel(expr):
    return infer_property_test_function(expr, properties.ROW_PANEL, isRowPanelTF)

def isRowPanelTF(node):
    if isinstance(node, ae.Scalar):
        return False
    else:
        size = node.size
        if size[0] < size[1]:
            return True
        else:
            return False

def isUnitary(node):
    if isinstance(node, ae.Symbol):
        return infer_property_symbol(node, properties.UNITARY, isUnitary)
    return False
        
def isNormal(node):
    if isinstance(node, ae.Symbol):
        return infer_property_symbol(node, properties.NORMAL, isNormal)
    return False
        
def isHermitian(node):
    if isinstance(node, ae.Symbol):
        return infer_property_symbol(node, properties.HERMITIAN, isHermitian)
    return False
        
def isPermutation(node):
    if isinstance(node, ae.Symbol):
        return infer_property_symbol(node, properties.PERMUTATION, isPermutation)
    if isinstance(node, ae.Times):
        return all(isPermutation(factor) for factor in node.operands)
    if isinstance(node, ae.Transpose):
        return isPermutation(node.operand)
    if isinstance(node, ae.Inverse):
        return isPermutation(node.operand)
    return False

def admitsFactorization(expr):
    return infer_property_test_function(expr, properties.ADMITS_FACTORIZATION, admitsFactorizationTF)

def admitsFactorizationTF(node):
    return not (isDiagonalB(node) or isTriangularB(node) or isOrthogonal(node) or isOrthogonalColumns(node) or isOrthogonalRows(node) or isPermutation(node))

def isScalar(expr):
    return infer_property_test_function(expr, properties.SCALAR, isScalarTF)

def isScalarTF(node):
    size = node.size
    if size[0] == 1 and size[1] == 1:
        return True
    else:
        return False

def isVector(expr):
    return infer_property_test_function(expr, properties.VECTOR, isVectorTF)

def isVectorTF(node):
    rows, columns = node.size
    if (rows == 1 and columns != 1) or (rows != 1 and columns == 1):
        return True
    else:
        return False

def isMatrix(expr):
    return infer_property_test_function(expr, properties.MATRIX, isMatrixTF)

def isMatrixTF(node):
    size = node.size
    if size[0] != 1 and size[1] != 1:
        return True
    else:
        return False

def isLowerTriangular(expr):
    return infer_property_test_function(expr, properties.LOWER_TRIANGULAR, isLowerTriangularB)

def isLowerTriangularB(node):
    lb, ub = node.bandwidth
    if ub <= 0:
        return True
    else:
        return False

def isUpperTriangular(expr):
    return infer_property_test_function(expr, properties.UPPER_TRIANGULAR, isUpperTriangularB)

def isUpperTriangularB(node):
    lb, ub = node.bandwidth
    if lb <= 0:
        return True
    else:
        return False

def isTriangular(expr):
    return infer_property_test_function(expr, properties.TRIANGULAR, isTriangularB)

def isTriangularB(node):
    lb, ub = node.bandwidth
    if lb <= 0 or ub <= 0:
        return True
    else:
        return False

def isDiagonal(expr):
    return infer_property_test_function(expr, properties.DIAGONAL, isDiagonalB)

def isDiagonalB(node):
    lb, ub = node.bandwidth
    if lb == 0 and ub == 0:
        return True
    else:
        return False

def isZero(expr):
    return infer_property_test_function(expr, properties.ZERO, isZeroB)

def isZeroB(node):
    # Careful. Unless the bandwidth of the Zero matrix is set to something
    # appropriate (e.g. (0, -1)), this is not a property that purely depends
    # on bandwidth.
    lb, ub = node.bandwidth
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
    properties.SYMMETRIC: isSymmetric,
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
    # properties.EXISTS_LU: isExistsLU,
    properties.ADMITS_FACTORIZATION: admitsFactorization,
    properties.FACTOR: isFactor,
    properties.CONSTANT: isConstant,
}

def infer_property(expr, prop):
    if not expr.is_constant:
        return False
    try:
        func = property_to_function[prop]
    except KeyError:
        raise NotImplementedError("No function to infer %s" % prop)
    return func(expr)

if __name__ == "__main__":
    pass