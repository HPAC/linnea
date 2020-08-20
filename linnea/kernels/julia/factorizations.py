from ...algebra.expression import Times, Plus, Transpose, \
                                  Scalar, Matrix

from ..utils.general import SizeArgument, \
                                       StrideArgument, \
                                       StorageFormatArgument, \
                                       InputOperand
                                       
from ..utils.factorizations import FactorizationKernel, OutputOperand

from ...utils import CodeTemplate, PropertyConstraint

from ...code_generation.memory.storage_format import StorageFormat

from ...algebra.properties import Property

from ...algebra import expression as ae

import matchpy


#####################
# Cholesky
#
# id 0
# Note: We assume that all symmetric matrices are stored as lower triangular
# matrices. Actually, that's not necessary in Julia, obtaining the other half is
# just more expensive. Test which conversion is more expensive.
# Actually, with storage format conversions, we don't need this assumption.

_A = matchpy.Wildcard.symbol("_A", symbol_type=ae.Matrix)
_L = matchpy.Wildcard.symbol("_L")
cf = lambda d: (d["N"]**3)/3

cholesky = FactorizationKernel(
    matchpy.Pattern(_A, PropertyConstraint("_A", {Property.SPSD})),
    [InputOperand(_A, StorageFormat.symmetric_triangular)],
    Times(_L, Transpose(_L)),
    [OutputOperand(_L, _A, ("N", "N"), [Property.LOWER_TRIANGULAR, Property.NON_SINGULAR], StorageFormat.lower_triangular)],
    cf,
    None,
    CodeTemplate("""LinearAlgebra.LAPACK.potrf!('L', $_A)"""),
    None,
    [SizeArgument("N", _A, "rows")],
    )


#####################
# LU with pivoting
#
# id 1

_A = matchpy.Wildcard.symbol("_A", symbol_type=ae.Matrix)
_P = matchpy.Wildcard.symbol("_P")
_L = matchpy.Wildcard.symbol("_L")
_U = matchpy.Wildcard.symbol("_U")
cf = lambda d: 2*(d["N"]**3)/3

plu = FactorizationKernel(
    matchpy.Pattern(_A, PropertyConstraint("_A", {Property.NON_SINGULAR})),
    [InputOperand(_A, StorageFormat.full)],
    Times(Transpose(_P), _L, _U),
    [OutputOperand(_L, _A, ("N", "N"), [Property.LOWER_TRIANGULAR, Property.UNIT_DIAGONAL, Property.NON_SINGULAR], StorageFormat.lower_triangular_udiag),
     OutputOperand(_U, _A, ("N", "N"), [Property.UPPER_TRIANGULAR, Property.NON_SINGULAR], StorageFormat.upper_triangular),
     OutputOperand(_P, None, ("N", "N"), [Property.PERMUTATION], StorageFormat.ipiv)
    ],
    cf,
    CodeTemplate(),
    CodeTemplate("($_A, $_P, info) = LinearAlgebra.LAPACK.getrf!($_A)"),
    CodeTemplate(),
    [SizeArgument("N", _A, "rows")],
    )


#####################
# QR square
#
# id 2

# TODO update indices in the generation of replacment rules

_A = matchpy.Wildcard.symbol("_A", symbol_type=ae.Matrix)
_Q = matchpy.Wildcard.symbol("_Q")
_R = matchpy.Wildcard.symbol("_R")

def cf_qr_square(d):
    cost = 2*d["N"]**2*(d["M"]-d["N"]/3) # for R
    cost += 4*(d["M"]**2*d["N"]-d["M"]*d["N"]**2+d["N"]**3/3)
    return cost

# TODO I think A does not have to be a full rank matrix. In that case, I don't get a full Q.

qr_square = FactorizationKernel(
    matchpy.Pattern(_A, PropertyConstraint("_A", {Property.FULL_RANK, Property.SQUARE})),
    [InputOperand(_A, StorageFormat.full)],
    Times(_Q, _R),
    [OutputOperand(_Q, _A, ("N", "N"), [Property.FULL_RANK, Property.SQUARE, Property.ORTHOGONAL], StorageFormat.QRfact_Q),
     OutputOperand(_R, _A, ("N", "N"), [Property.UPPER_TRIANGULAR, Property.SQUARE], StorageFormat.QRfact_R)
    ],
    cf_qr_square,
    CodeTemplate(),
    CodeTemplate("$_A = qr!($_A)"),
    CodeTemplate(),
    [SizeArgument("M", _A, "rows"),
     SizeArgument("N", _A, "columns")],
    )


#####################
# QR column panel
#
# id 3

_A = matchpy.Wildcard.symbol("_A", symbol_type=ae.Matrix)
_Q = matchpy.Wildcard.symbol("_Q")
_R = matchpy.Wildcard.symbol("_R")

def cf_qr_column(d):
    cost = 2*d["N"]**2*(d["M"]-d["N"]/3) # for R
    cost += 2*d["N"]**2*(d["M"]-d["N"]/3) # for Q (because it has size m x n)
    return cost

# TODO I think A does not have to be a full rank matrix. In that case, I don't get a full Q.

qr_column = FactorizationKernel(
    matchpy.Pattern(_A, PropertyConstraint("_A", {Property.FULL_RANK, Property.COLUMN_PANEL})),
    [InputOperand(_A, StorageFormat.full)],
    Times(_Q, _R),
    [OutputOperand(_Q, _A, ("M", "N"), [Property.FULL_RANK, Property.COLUMN_PANEL, Property.ORTHOGONAL_COLUMNS], StorageFormat.QRfact_Q),
     OutputOperand(_R, _A, ("N", "N"), [Property.UPPER_TRIANGULAR, Property.SQUARE], StorageFormat.QRfact_R)
    ],
    cf_qr_column,
    None,
    CodeTemplate("$_A = qr!($_A)"),
    None,
    [SizeArgument("M", _A, "rows"),
     SizeArgument("N", _A, "columns")],
    )

#####################
# Eigendecomposition
#
# id 4

_A = matchpy.Wildcard.symbol("_A", symbol_type=ae.Matrix)
_Z = matchpy.Wildcard.symbol("_Z")
_W = matchpy.Wildcard.symbol("_W")

def cf_eigen(d):
    cost = 4/3*d["N"]**3 # reduction to tridiagional form
    cost += d["N"]**2 # solve
    cost += 2*d["N"]**3 # backtransformation
    return cost

eigendecomposition = FactorizationKernel(
    matchpy.Pattern(_A, PropertyConstraint("_A", {Property.SYMMETRIC})),
    [InputOperand(_A, StorageFormat.symmetric_triangular)],
    Times(_Z, _W, Transpose(_Z)),
    [OutputOperand(_Z, _A, ("N", "N"), [Property.ORTHOGONAL], StorageFormat.full),
     OutputOperand(_W, None, ("N", "N"), [Property.DIAGONAL, Property.SQUARE], StorageFormat.diagonal_vector)
    ],
    cf_eigen,
    None,
    CodeTemplate("$_W, $_A = LinearAlgebra.LAPACK.syev!('V', $uplo, $_A)"),
    None,
    [SizeArgument("N", _A, "rows"),
     StorageFormatArgument("uplo", _A, StorageFormat.symmetric_lower_triangular, ["L", "U"]),]
    )


#####################
# Singular value decomposition
#
# id 5

# TODO check number of FLOPs of economy SVD

# SVD column panel

_A = matchpy.Wildcard.symbol("_A", symbol_type=ae.Matrix)
_U = matchpy.Wildcard.symbol("_U")
_S = matchpy.Wildcard.symbol("_S")
_V = matchpy.Wildcard.symbol("_V")
cf = lambda d: 14*d["M"]*d["N"]**2 + 8*d["N"]**3

singular_value_cp = FactorizationKernel(
    matchpy.Pattern(_A, PropertyConstraint("_A", {Property.COLUMN_PANEL})),
    [InputOperand(_A, StorageFormat.full)],
    Times(_U, _S, _V),
    [OutputOperand(_U, _A, ("M", "N"), [Property.ORTHOGONAL_COLUMNS], StorageFormat.full),
     OutputOperand(_S, None, ("N", "N"), [Property.DIAGONAL], StorageFormat.diagonal_vector),
     OutputOperand(_V, None, ("N", "N"), [Property.ORTHOGONAL], StorageFormat.full)
    ],
    cf,
    None,
    CodeTemplate("(_, $_S, $_V) = LAPACK.gesvd!('O', 'S', $_A)"),
    None,
    [SizeArgument("M", _A, "rows"),
     SizeArgument("N", _A, "columns")],
    )

# SVD square
# id 6

_A = matchpy.Wildcard.symbol("_A", symbol_type=ae.Matrix)
_U = matchpy.Wildcard.symbol("_U")
_S = matchpy.Wildcard.symbol("_S")
_V = matchpy.Wildcard.symbol("_V")
cf = lambda d: 22*d["N"]**3

singular_value_sq = FactorizationKernel(
    matchpy.Pattern(_A, PropertyConstraint("_A", {Property.SQUARE})),
    [InputOperand(_A, StorageFormat.full)],
    Times(_U, _S, _V),
    [OutputOperand(_U, _A, ("N", "N"), [Property.ORTHOGONAL], StorageFormat.full),
     OutputOperand(_S, None, ("N", "N"), [Property.DIAGONAL], StorageFormat.diagonal_vector),
     OutputOperand(_V, None, ("N", "N"), [Property.ORTHOGONAL], StorageFormat.full)
    ],
    cf,
    None,
    CodeTemplate("(_, $_S, $_V) = LAPACK.gesvd!('O', 'S', $_A)"),
    None,
    [SizeArgument("N", _A, "rows")],
    )

# SVD row panel
# id 7

_A = matchpy.Wildcard.symbol("_A", symbol_type=ae.Matrix)
_U = matchpy.Wildcard.symbol("_U")
_S = matchpy.Wildcard.symbol("_S")
_V = matchpy.Wildcard.symbol("_V")
cf = lambda d: 14*d["M"]*d["N"]**2 + 8*d["N"]**3

singular_value_rp = FactorizationKernel(
    matchpy.Pattern(_A, PropertyConstraint("_A", {Property.ROW_PANEL})),
    [InputOperand(_A, StorageFormat.full)],
    Times(_U, _S, _V),
    [OutputOperand(_U, None, ("M", "M"), [Property.ORTHOGONAL], StorageFormat.full),
     OutputOperand(_S, None, ("M", "M"), [Property.DIAGONAL], StorageFormat.diagonal_vector),
     OutputOperand(_V, _A, ("M", "N"), [Property.ORTHOGONAL_ROWS], StorageFormat.full)
    ],
    cf,
    None,
    CodeTemplate("($_U, $_S, _) = LAPACK.gesvd!('S', 'O', $_A)"),
    None,
    [SizeArgument("M", _A, "rows"),
     SizeArgument("N", _A, "columns")],
    )