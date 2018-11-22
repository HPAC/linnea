import enum

from .. import utils

@utils.PartiallyOrderedEnum
class Property(enum.Enum):
    INPUT = "Input"
    OUTPUT = "Output"
    AUXILIARY = "Auxiliary"
    #
    SCALAR = "Scalar"
    VECTOR = "Vector"
    MATRIX = "Matrix"
    #
    SQUARE       = "Square"
    ROW_PANEL    = "RowPanel"
    COLUMN_PANEL = "ColumnPanel"
    #
    ZERO = "Zero"
    IDENTITY = "Identity"
    CONSTANT = "Constant"
    DIAGONAL = "Diagonal"
    TRIANGULAR = "Triangular"
    LOWER_TRIANGULAR = "LowerTriangular"
    UPPER_TRIANGULAR = "UpperTriangular"
    UNIT_DIAGONAL = "UnitDiagonal"
    #
    POSITIVE = "Positive" # only for scalars
    #
    SYMMETRIC = "Symmetric"
    SPSD = "SPSD"
    SPD = "SPD"
    #
    NON_SINGULAR = "Non-singular"
    ORTHOGONAL   = "Orthogonal"
    FULL_RANK    = "FullRank"
    #
    UNITARY      = "Unitary"
    NORMAL       = "Normal"
    HERMITIAN    = "Hermitian"
    #
    ORTHOGONAL_COLUMNS   = "OrthogonalColumns" # for QR
    ORTHOGONAL_ROWS   = "OrthogonalRows" # for SVD
    #
    PERMUTATION = "Permutation"
    #
    ADMITS_FACTORIZATION = "AdmitsFactorization"
    #
    FACTOR = "Factor"

    __ordering__ = {
        (ZERO, CONSTANT),
        (IDENTITY, CONSTANT),
        (IDENTITY, DIAGONAL),
        (IDENTITY, FULL_RANK),
        (DIAGONAL, LOWER_TRIANGULAR),
        (DIAGONAL, UPPER_TRIANGULAR),
        (LOWER_TRIANGULAR, TRIANGULAR),
        (UPPER_TRIANGULAR, TRIANGULAR),
        (SPD, SPSD),
        (SPSD, SYMMETRIC),
        (SYMMETRIC, SQUARE),
        (SPD, NON_SINGULAR),
        (PERMUTATION, ORTHOGONAL),
        (ORTHOGONAL, SQUARE),
        (ORTHOGONAL, ORTHOGONAL_ROWS), # this can be used to simplify the simplify() function
        (ORTHOGONAL, ORTHOGONAL_COLUMNS), # this can be used to simplify the simplify() function
        (ORTHOGONAL_ROWS, FULL_RANK),
        (ORTHOGONAL_COLUMNS, FULL_RANK),
        (ORTHOGONAL, NON_SINGULAR),
        (NON_SINGULAR, FULL_RANK)
    }


# Property.to_dot_file(Property)


# Add dictionary implies (prop) -> props
#   e.g., SPD -> Symmetric
#         Diagonal -> Triangular
# RHS has to be list of properties

# TODO what about "negative" implications?
#   Square -> not RowPanel

# TODO what about coobinations of properties?
#   lower and upper triangular -> diagonal

implications = dict()
for p1, p2 in Property.__transitive_closure__:
    implications.setdefault(Property(p1), set()).add(Property(p2))


# mutually_exclusive = [
#     set((Property.INPUT, Property.OUTPUT, Property.AUXILIARY)),
#     set((Property.MATRIX, Property.VECTOR. Property.SCALAR)),
#     set((Property.SQUARE, Property.ROW_PANEL, Property.COLUMN_PANEL)),
#     set((Property.IDENTITY, Property.ZERO)),
# ]

negative_implications = {
    Property.INPUT : [Property.OUTPUT, Property.AUXILIARY],
    Property.OUTPUT : [Property.INPUT, Property.AUXILIARY],
    Property.AUXILIARY : [Property.INPUT, Property.OUTPUT],
    #
    Property.SCALAR : [Property.VECTOR, Property.MATRIX],
    Property.VECTOR : [Property.SCALAR, Property.MATRIX],
    Property.MATRIX : [Property.SCALAR, Property.VECTOR],
    #
    Property.SQUARE       : [Property.ROW_PANEL, Property.COLUMN_PANEL],
    Property.ROW_PANEL    : [Property.SQUARE],
    Property.COLUMN_PANEL : [Property.SQUARE],
    #
    Property.ZERO : [Property.IDENTITY, Property.NON_SINGULAR, Property.FULL_RANK],
    Property.IDENTITY : [Property.ROW_PANEL, Property.COLUMN_PANEL, Property.ZERO],
    Property.CONSTANT : [],
    Property.DIAGONAL : [],
    Property.TRIANGULAR : [], # is triangular a contradiction to symmetric? no if matrix is diagonal
    Property.LOWER_TRIANGULAR : [],
    Property.UPPER_TRIANGULAR : [],
    Property.UNIT_DIAGONAL : [], # TODO missing
    #
    Property.POSITIVE : [],
    Property.SYMMETRIC : [Property.ROW_PANEL, Property.COLUMN_PANEL],
    Property.SPSD : [Property.ROW_PANEL, Property.COLUMN_PANEL],
    Property.SPD : [Property.ROW_PANEL, Property.COLUMN_PANEL],
    #
    Property.NON_SINGULAR : [Property.ROW_PANEL, Property.COLUMN_PANEL],
    Property.ORTHOGONAL   : [Property.ROW_PANEL, Property.COLUMN_PANEL],
    Property.FULL_RANK    : [], # TODO missing
    # EXISTS_LU    : [], # TODO missing
    #
    Property.UNITARY      : [Property.ROW_PANEL, Property.COLUMN_PANEL],
    Property.NORMAL       : [Property.ROW_PANEL, Property.COLUMN_PANEL],
    Property.HERMITIAN    : [Property.ROW_PANEL, Property.COLUMN_PANEL],
    #
    Property.ORTHOGONAL_COLUMNS : [],
    Property.ORTHOGONAL_ROWS : [],
    #
    Property.PERMUTATION : [Property.ROW_PANEL, Property.COLUMN_PANEL],
    #
    Property.ADMITS_FACTORIZATION : [],
}

if __name__ == "__main__":
    pass