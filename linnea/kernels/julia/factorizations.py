from ...algebra.expression import Times, Plus, Transpose, \
                                  Scalar, Matrix

# from patternmatcher.functional import Pattern

from ..utils.general import SizeArgument, \
                                       StrideArgument, \
                                       StorageFormatArgument, \
                                       InputOperand
                                       
from ..utils.factorizations import FactorizationKernel, OutputOperand

from ...utils import CodeTemplate

from ...code_generation.memory.storage_format import StorageFormat

from ...algebra.properties import Property as properties

from ...algebra import expression as ae

import matchpy


#####################
# Cholesky
#
# Note: We assume that all symmetric matrices are stored as lower triangular
# matrices. Actually, that's not necessary in Julia, obtaining the other half is
# just more expensive. Test which conversion is more expensive.
# Actually, with storage format conversions, we don't need this assumption.

_A = matchpy.Wildcard.symbol("_A", symbol_type=ae.Matrix)
_L = matchpy.Wildcard.symbol("_L")
cf = lambda d: (d["N"]**3)/3

cholesky = FactorizationKernel(
    matchpy.Pattern(_A,
                    matchpy.CustomConstraint(lambda _A: _A.has_property(properties.SPD))
                    ),
    [InputOperand(_A, StorageFormat.symmetric_triangular)],
    Times(_L, Transpose(_L)),
    [OutputOperand(_L, _A, ("N", "N"), [properties.LOWER_TRIANGULAR, properties.NON_SINGULAR], StorageFormat.lower_triangular)],
    cf,
    None,
    CodeTemplate("""Base.LinAlg.LAPACK.potrf!('L', $_A)"""),
    # CodeTemplate("""$_A = cholfact!($_A, :L)"""),
    None,
    [SizeArgument("N", _A, "rows")],
    )

#####################
# LU with pivoting
#

_A = matchpy.Wildcard.symbol("_A", symbol_type=ae.Matrix)
_P = matchpy.Wildcard.symbol("_P")
_L = matchpy.Wildcard.symbol("_L")
_U = matchpy.Wildcard.symbol("_U")
cf = lambda d: 2*(d["N"]**2)/3

plu = FactorizationKernel(
    matchpy.Pattern(_A,
                    matchpy.CustomConstraint(
                        lambda _A: _A.has_property(properties.NON_SINGULAR)
                    )
    ),
    [InputOperand(_A, StorageFormat.full)],
    Times(Transpose(_P), _L, _U),
    [OutputOperand(_L, _A, ("N", "N"), [properties.LOWER_TRIANGULAR, properties.UNIT_DIAGONAL, properties.NON_SINGULAR], StorageFormat.LUfact_L),
     OutputOperand(_U, _A, ("N", "N"), [properties.UPPER_TRIANGULAR, properties.NON_SINGULAR], StorageFormat.LUfact_U),
     OutputOperand(_P, _A, ("N", "N"), [properties.PERMUTATION], StorageFormat.LUfact_P)
    ],
    cf,
    CodeTemplate(),
    CodeTemplate("$_A = lufact!($_A)"),
    CodeTemplate(),
    [SizeArgument("N", _A, "rows")],
    )


#####################
# QR square
#

_A = matchpy.Wildcard.symbol("_A", symbol_type=ae.Matrix)
_Q = matchpy.Wildcard.symbol("_Q")
_R = matchpy.Wildcard.symbol("_R")

def cf_qr_square(d):
    cost = 2*d["N"]**2*(d["M"]-d["N"]/3) # for R
    cost += 4*(d["M"]**2*d["N"]-d["M"]*d["N"]**2+d["N"]**3/3)
    return cost

# TODO I think A does not have to be a full rank matrix. In that case, I don't get a full Q.

qr_square = FactorizationKernel(
    matchpy.Pattern(_A,
                    matchpy.CustomConstraint(
                        lambda _A: _A.has_property(properties.FULL_RANK) and \
                              _A.has_property(properties.SQUARE)
                    )
    ),
    [InputOperand(_A, StorageFormat.full)],
    Times(_Q, _R),
    [OutputOperand(_Q, _A, ("N", "N"), [properties.FULL_RANK, properties.SQUARE, properties.ORTHOGONAL], StorageFormat.QRfact_Q),
     OutputOperand(_R, _A, ("N", "N"), [properties.UPPER_TRIANGULAR, properties.SQUARE], StorageFormat.QRfact_R)
    ],
    cf_qr_square,
    CodeTemplate(),
    CodeTemplate("$_A = qrfact!($_A)"),
    CodeTemplate(),
    [SizeArgument("M", _A, "rows"),
     SizeArgument("N", _A, "columns")],
    )


#####################
# QR column panel
#

_A = matchpy.Wildcard.symbol("_A", symbol_type=ae.Matrix)
_Q = matchpy.Wildcard.symbol("_Q")
_R = matchpy.Wildcard.symbol("_R")

def cf_qr_column(d):
    cost = 2*d["N"]**2*(d["M"]-d["N"]/3) # for R
    cost += 2*d["N"]**2*(d["M"]-d["N"]/3) # for Q (because it has size m x n)
    return cost

# TODO I think A does not have to be a full rank matrix. In that case, I don't get a full Q.

qr_column = FactorizationKernel(
    matchpy.Pattern(_A,
                    matchpy.CustomConstraint(
                        lambda _A: _A.has_property(properties.FULL_RANK) and \
                          _A.has_property(properties.COLUMN_PANEL)
                    )
    ),
    [InputOperand(_A, StorageFormat.full)],
    Times(_Q, _R),
    [OutputOperand(_Q, _A, ("M", "N"), [properties.FULL_RANK, properties.COLUMN_PANEL, properties.ORTHOGONAL_COLUMNS], StorageFormat.QRfact_Q),
     OutputOperand(_R, _A, ("N", "N"), [properties.UPPER_TRIANGULAR, properties.SQUARE], StorageFormat.QRfact_R)
    ],
    cf_qr_column,
    None,
    CodeTemplate("$_A = qrfact!($_A)"),
    None,
    [SizeArgument("M", _A, "rows"),
     SizeArgument("N", _A, "columns")],
    )

#####################
# Eigendecomposition
#
# syev
# probably no problem

_A = matchpy.Wildcard.symbol("_A", symbol_type=ae.Matrix)
_Z = matchpy.Wildcard.symbol("_Z")
_W = matchpy.Wildcard.symbol("_W")

def cf_eigen(d):
    cost = 4/3*d["N"]**3 # reduction to tridiagional form
    cost += d["N"]**2 # solve
    cost += 2*d["N"]**3 # backtransformation
    return cost

eigendecomposition = FactorizationKernel(
    matchpy.Pattern(_A,
                    matchpy.CustomConstraint(
                        lambda _A: _A.has_property(properties.SYMMETRIC)
                    )
    ),
    [InputOperand(_A, StorageFormat.symmetric_triangular)],
    Times(_Z, _W, Transpose(_Z)),
    [OutputOperand(_Z, _A, ("N", "N"), [properties.ORTHOGONAL], StorageFormat.full),
     OutputOperand(_W, None, ("N", "N"), [properties.DIAGONAL, properties.SQUARE], StorageFormat.diagonal_vector)
    ],
    cf_eigen,
    None,
    CodeTemplate("$_W, $_A = Base.LinAlg.LAPACK.syev!('V', $uplo, $_A)"),
    # CodeTemplate("$_A = eigfact!($_A)"),
    None,
    [SizeArgument("N", _A, "rows"),
     StorageFormatArgument("uplo", _A, StorageFormat.symmetric_lower_triangular, ["L", "U"]),]
    )


#####################
# Singular value decomposition
#
# gesvd
# probably no problem

_A = matchpy.Wildcard.symbol("_A", symbol_type=ae.Matrix)
_U = matchpy.Wildcard.symbol("_U")
_S = matchpy.Wildcard.symbol("_S")
_V = matchpy.Wildcard.symbol("_V")
cf = lambda d: 14*d["M"]*d["N"]**2 + 8*d["N"]**3

singular_value = FactorizationKernel(
    matchpy.Pattern(_A),
    [InputOperand(_A, StorageFormat.full)],
    Times(_U, _S, _V),
    [OutputOperand(_U, _A, ("M", "M"), [properties.ORTHOGONAL], StorageFormat.svdfact_U),
     OutputOperand(_S, _A, ("M", "N"), [properties.DIAGONAL], StorageFormat.svdfact_S), # for economy SVD, sizes do not depend on rank. Remember to simplify storage format transformation.
     OutputOperand(_V, _A, ("N", "N"), [properties.ORTHOGONAL], StorageFormat.svdfact_V) # this changes to OrthogonalRows for economy SVD
    ],
    cf,
    None,
    CodeTemplate("$_A = svdfact!($_A, thin=false)"),
    None,
    [SizeArgument("M", _A, "rows"),
     SizeArgument("N", _A, "columns")],
    )

