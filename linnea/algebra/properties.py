import enum

from .. import utils

class PropertyError(Exception):
    pass

@utils.PartiallyOrderedEnum
class Property(enum.Enum):
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



implications = dict()
for p1, p2 in Property.__transitive_closure__:
    implications.setdefault(Property(p1), set()).add(Property(p2))


binary_implications = {
    # This is already covered by how bandwidth is set.
    # frozenset({Property.LOWER_TRIANGULAR, Property.UPPER_TRIANGULAR}): Property.DIAGONAL,
    frozenset({Property.SQUARE, Property.DIAGONAL}): Property.SYMMETRIC, # ESSENTIAL
    # Orthogonal rows and columns is only possible if the matrix is square. Redundant?
    # frozenset({Property.ORTHOGONAL_ROWS, Property.ORTHOGONAL_COLUMNS}): Property.ORTHOGONAL,
    # frozenset({Property.ORTHOGONAL_ROWS, Property.SQUARE}): Property.ORTHOGONAL,
    # frozenset({Property.ORTHOGONAL_COLUMNS, Property.SQUARE}): Property.ORTHOGONAL,
    # frozenset({Property.SPSD, Property.NON_SINGULAR}): Property.SPD,
    # TODO problem: key is the same here. Right now, those four are covered by construction of IdentityMatrix.
    # frozenset({Property.IDENTITY, Property.SQUARE}): Property.SPD, # ESSENTIAL
    # frozenset({Property.IDENTITY, Property.SQUARE}): Property.PERMUTATION, # ESSENTIAL
    # frozenset({Property.IDENTITY, Property.COLUMN_PANEL}): Property.ORTHOGONAL_COLUMNS,
    # frozenset({Property.IDENTITY, Property.ROW_PANEL}): Property.ORTHOGONAL_ROWS,
    frozenset({Property.FULL_RANK, Property.SQUARE}): Property.NON_SINGULAR, # ESSENTIAL
}

binary_implications_backwards = dict()
for s, p in binary_implications.items():
    binary_implications_backwards.setdefault(p, []).append(s)

# for k, v in binary_implications_backwards.items():
#     print(k, v)
# print(binary_implications_backwards)

negative_implications = {
    Property.AUXILIARY: set(),
    #
    Property.SCALAR: {Property.VECTOR, Property.MATRIX},
    Property.VECTOR: {Property.SCALAR, Property.MATRIX,
                      Property.SYMMETRIC, Property.SPSD, Property.SPD},
    Property.MATRIX: {Property.SCALAR, Property.VECTOR},
    #
    Property.SQUARE      : {Property.ROW_PANEL, Property.COLUMN_PANEL},
    Property.ROW_PANEL   : {Property.SQUARE, Property.COLUMN_PANEL,
                            Property.ORTHOGONAL_COLUMNS,
                            Property.SYMMETRIC, Property.SPSD, Property.SPD},
    Property.COLUMN_PANEL: {Property.SQUARE, Property.ROW_PANEL,
                            Property.ORTHOGONAL_ROWS,
                            Property.SYMMETRIC, Property.SPSD, Property.SPD},
    #
    Property.ZERO            : {Property.IDENTITY, Property.NON_SINGULAR,
                                Property.FULL_RANK},
    Property.IDENTITY        : {Property.ZERO},
    Property.CONSTANT        : set(),
    #
    Property.DIAGONAL        : set(),
    Property.TRIANGULAR      : set(),
    Property.LOWER_TRIANGULAR: set(),
    Property.UPPER_TRIANGULAR: set(),
    Property.UNIT_DIAGONAL   : set(), # TODO missing
    #
    Property.POSITIVE : {Property.VECTOR, Property.MATRIX},
    Property.SYMMETRIC: {Property.ROW_PANEL, Property.COLUMN_PANEL,
                         Property.VECTOR},
    Property.SPSD     : {Property.ROW_PANEL, Property.COLUMN_PANEL,
                         Property.VECTOR},
    Property.SPD      : {Property.ROW_PANEL, Property.COLUMN_PANEL,
                         Property.VECTOR},
    #
    Property.NON_SINGULAR: {Property.ROW_PANEL, Property.COLUMN_PANEL,
                            Property.ZERO},
    Property.FULL_RANK   : {Property.ZERO},
    #
    Property.ORTHOGONAL        : {Property.ROW_PANEL, Property.COLUMN_PANEL,
                                  Property.ZERO},
    Property.ORTHOGONAL_COLUMNS: set(),
    Property.ORTHOGONAL_ROWS   : set(),
    Property.PERMUTATION       : {Property.ROW_PANEL, Property.COLUMN_PANEL,
                                  Property.ZERO},
    #
    Property.UNITARY     : {Property.ROW_PANEL, Property.COLUMN_PANEL},
    Property.NORMAL      : {Property.ROW_PANEL, Property.COLUMN_PANEL},
    Property.HERMITIAN   : {Property.ROW_PANEL, Property.COLUMN_PANEL},
    #
    Property.ADMITS_FACTORIZATION: {Property.DIAGONAL, Property.TRIANGULAR,
                                    Property.LOWER_TRIANGULAR,
                                    Property.UPPER_TRIANGULAR,
                                    Property.ORTHOGONAL,
                                    Property.ORTHOGONAL_COLUMNS,
                                    Property.ORTHOGONAL_ROWS,
                                    Property.PERMUTATION, Property.IDENTITY,
                                    Property.ZERO},
}

if __name__ == "__main__":
    pass