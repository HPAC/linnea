import enum

from .. import utils

@utils.PartiallyOrderedEnum
class Property(enum.Enum):
    # Input = enum.auto()
    # Output = enum.auto()
    # Auxiliary = enum.auto()
    # #
    # Scalar = enum.auto()
    # Vector = enum.auto()
    # Matrix = enum.auto()
    # #
    # Square = enum.auto()
    # RowPanel = enum.auto()
    # ColumnPanel = enum.auto()
    # #
    # Zero = enum.auto()
    # Identity = enum.auto()
    # Diagonal = enum.auto()
    # Triangular = enum.auto()
    # LowerTriangular = enum.auto()
    # UpperTriangular = enum.auto()
    # UnitDiagonal = enum.auto()
    # #
    # Symmetric = enum.auto()
    # SPD = enum.auto()
    # #
    # NonSingular = enum.auto()
    # Orthogonal = enum.auto()
    # FullRank = enum.auto()
    # #
    # Unitary = enum.auto()
    # Normal = enum.auto()
    # Hermitian = enum.auto()
    # #
    # OrthogonalColumns = enum.auto() # for QR
    # #
    # Permutation = enum.auto()
    # #
    # AdmitsFactorization = enum.auto()

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
    SYMMETRIC = "Symmetric"
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
        # (IDENTITY, SYMMETRIC),
        (IDENTITY, FULL_RANK),
        # (IDENTITY, NON_SINGULAR),
        (DIAGONAL, LOWER_TRIANGULAR),
        (DIAGONAL, UPPER_TRIANGULAR),
        (LOWER_TRIANGULAR, TRIANGULAR),
        (UPPER_TRIANGULAR, TRIANGULAR),
        (SPD, SYMMETRIC),
        (SYMMETRIC, SQUARE),
        (SPD, NON_SINGULAR),
        # (SPD, FULL_RANK),
        (PERMUTATION, ORTHOGONAL),
        (ORTHOGONAL, SQUARE),
        (ORTHOGONAL, ORTHOGONAL_ROWS), # this can be used to simplify the simplify() function
        (ORTHOGONAL, ORTHOGONAL_COLUMNS), # this can be used to simplify the simplify() function
        (ORTHOGONAL_ROWS, FULL_RANK),
        (ORTHOGONAL_COLUMNS, FULL_RANK),
        (ORTHOGONAL, NON_SINGULAR),
        # (ORTHOGONAL_ROWS, NON_SINGULAR),
        # (ORTHOGONAL_COLUMNS, NON_SINGULAR),
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

# print(implications)

# implications = {
#     INPUT : [],
#     OUTPUT : [],
#     AUXILIARY : [],
#     #
#     SCALAR : [],
#     VECTOR : [],
#     MATRIX : [],
#     #
#     SQUARE       : [],
#     ROW_PANEL    : [],
#     COLUMN_PANEL : [],
#     #
#     ZERO : [],
#     IDENTITY : [SQUARE, NON_SINGULAR, FULL_RANK, DIAGONAL],
#     DIAGONAL : [TRIANGULAR, LOWER_TRIANGULAR, UPPER_TRIANGULAR],
#     TRIANGULAR : [],
#     LOWER_TRIANGULAR : [TRIANGULAR],
#     UPPER_TRIANGULAR : [TRIANGULAR],
#     UNIT_DIAGONAL : [], # TODO missing
#     #
#     SYMMETRIC : [SQUARE],
#     SPD : [SQUARE, SYMMETRIC, NON_SINGULAR, FULL_RANK],
#     #
#     NON_SINGULAR : [FULL_RANK],
#     ORTHOGONAL   : [SQUARE, NON_SINGULAR, FULL_RANK],
#     FULL_RANK    : [], # TODO missing
#     # EXISTS_LU    : [], # TODO missing
#     #
#     UNITARY      : [SQUARE, NORMAL, NON_SINGULAR, FULL_RANK],
#     NORMAL       : [SQUARE],
#     HERMITIAN    : [SQUARE, NORMAL],
#     #
#     ORTHOGONAL_COLUMNS : [], # TODO missing
#     #
#     PERMUTATION : [SQUARE, NON_SINGULAR, ORTHOGONAL, FULL_RANK],
# }

# # TODO this looks way too simple

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
    Property.SYMMETRIC : [Property.ROW_PANEL, Property.COLUMN_PANEL],
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