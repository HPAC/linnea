# from patternmatcher.expression import Times, Plus, Transpose, \
#                                       Scalar, Matrix, \
#                                       WildcardSymbol

# from patternmatcher.functional import Pattern

# from clak.kernels.utils.general import CodeTemplate, \
#                                SizeArgument, StrideArgument
# from clak.kernels.utils.factorizations import FactorizationKernel, OutputOperand

# import patternmatcher.properties as properties

import matchpy

n = 10
m = 20


#####################
# Cholesky
#
# potrf
# Note: We assume that all symmetric matrices are stored as lower triangular
# matrices.

_A = WildcardSymbol("_A")
_L = WildcardSymbol("_L")
L = Matrix("L", (n, n))
cf = lambda d: d["N"]**3/3

cholesky = FactorizationKernel(
    Pattern(
        _A,
        lambda d: d["_A"].has_property(properties.SPD)
        ),
    Times(_L, Transpose(_L)),
    [OutputOperand(_L, _A, ("N", "N"), [properties.LOWER_TRIANGULAR])],
    cf,
    None,
    CodeTemplate("""${type_prefix}potrf("L", $N, $_A, $ldA, &info); // cholesky"""),
    None,# CodeTemplate("after\n")
    [SizeArgument("N", _A, "rows"),
     StrideArgument("ldA", _A, "rows")],
    )


#####################
# LDL
#
# sytrf
# There is no function to recunstruct L. Only solve with LDL and invert with LDL

_A = WildcardSymbol("_A")
_L = WildcardSymbol("_L")
_D = WildcardSymbol("_D")
L = Matrix("L", (n, n))
D = Matrix("D", (n, n))
cf = lambda d: 0 # TODO missing n^3/3 (Golub, van Loan)

# TODO overwriting missing

ldl = FactorizationKernel(
    Pattern(
        _A,
        lambda d: d["_A"].has_property(properties.SYMMETRIC) and \
                  not d["_A"].has_property(properties.DIAGONAL) and \
                  not d["_A"].has_property(properties.SPD)
        ),
    Times(_L, _D, Transpose(_L)),
    [OutputOperand(_L, None, ("N", "N"), [properties.LOWER_TRIANGULAR, properties.FULL_RANK, properties.NON_SINGULAR, properties.SQUARE]),
     OutputOperand(_D, None, ("N", "N"), [properties.DIAGONAL, properties.SQUARE])
    ],
    cf,
    CodeTemplate(),
    CodeTemplate("""${type_prefix}sytrf("L", $N, $_A, $ldA, $IPIV, $work, $lwork, &info); // LDL"""),
    CodeTemplate(),
    [SizeArgument("N", _A, "rows"),
     StrideArgument("ldA", _A, "rows")],
    )


#####################
# LU with pivoting
#
# getrf
# no problem

# TODO How to handle "multiplication" with pivoting matrix? Use (s/d)laswp
# What about (s/d)lapmt?

_A = WildcardSymbol("_A")
_P = WildcardSymbol("_P")
_L = WildcardSymbol("_L")
_U = WildcardSymbol("_U")
P = Matrix("P", (n, n))
L = Matrix("L", (n, n))
U = Matrix("U", (n, n))
cf = lambda d: 2*d["N"]**2/3

plu = FactorizationKernel(
    Pattern(
        _A,
        lambda d: d["_A"].has_property(properties.NON_SINGULAR) and \
                  not d["_A"].has_property(properties.DIAGONAL) and \
                  not d["_A"].has_property(properties.TRIANGULAR) and \
                  not d["_A"].has_property(properties.ORTHOGONAL) and \
                  not d["_A"].has_property(properties.SYMMETRIC)
        ),
    Times(_P, _L, _U),
    [OutputOperand(_L, _A, ("N", "N"), [properties.LOWER_TRIANGULAR, properties.UNIT_DIAGONAL]),
     OutputOperand(_U, _A, ("N", "N"), [properties.UPPER_TRIANGULAR]),
     OutputOperand(_P, None, ("N", "N"), [properties.PERMUTATION])
    ],
    cf,
    CodeTemplate(),
    CodeTemplate("${type_prefix}getrf($M, $N, $_A, $ldA, $_P, &info); // PLU"),
    CodeTemplate(),
    [SizeArgument("N", _A, "rows"),
     StrideArgument("ldA", _A, "rows")],
    )


#####################
# QR square
#
# geqrf
# Q is not generated explicitly, but can be reconstructed using orgqr
# Multiply by Q without constructing it: ormqr
# When trying to multiply with Q and ormqr can not be used, during code generation,
# test corresponding memory property and generate code to reconstruct Q. This
# could/has to be? tested before every operation. However, how to take advantage
# of ormqr? Problem when generating operation to reconstruct during code
# generation is that this changes the cost.
# Alternatively, implement that as real properties?

_A = WildcardSymbol("_A")
_Q = WildcardSymbol("_Q")
_R = WildcardSymbol("_R")
Q = Matrix("Q", (n, n))
R = Matrix("R", (n, n))

def cf_qr_square(d):
    cost = 2*d["N"]**2*(d["M"]-d["N"]/3) # for R
    cost += 4*(d["M"]**2*d["N"]-d["M"]*d["N"]**2+d["N"]**3/3)
    return cost

# TODO I think A does not have to be a full rank matrix. In that case, I don't get a full Q.

qr_square = FactorizationKernel(
    Pattern(
        _A,
        lambda d: d["_A"].has_property(properties.FULL_RANK) and \
                  d["_A"].has_property(properties.SQUARE) \
                  and not d["_A"].has_property(properties.TRIANGULAR) \
                  and not d["_A"].has_property(properties.DIAGONAL) \
                  and not d["_A"].has_property(properties.ORTHOGONAL) \
                  and not d["_A"].has_property(properties.ORTHOGONAL_COLUMNS)
        ),
    Times(_Q, _R),
    [OutputOperand(_Q, None, ("N", "N"), [properties.FULL_RANK, properties.SQUARE, properties.ORTHOGONAL]),
     OutputOperand(_R, _A, ("N", "N"), [properties.UPPER_TRIANGULAR, properties.SQUARE])
    ],
    cf_qr_square,
    CodeTemplate(),
    CodeTemplate("${type_prefix}geqrf($M, $N, $_A, $ldA, $tau, $work, $lwork, &info); // QR square"),
    CodeTemplate(),
    [SizeArgument("M", _A, "rows"),
     SizeArgument("N", _A, "columns"),
     StrideArgument("ldA", _A, "rows")],
    )


#####################
# QR column panel
#
# geqrf
# Q is not generated explicitly, but can be reconstructed using orgqr
# Multiply by Q without constructing it: ormqr

_A = WildcardSymbol("_A")
_Q = WildcardSymbol("_Q")
_R = WildcardSymbol("_R")
Q = Matrix("Q", (m, n))
R = Matrix("R", (n, n))

def cf_qr_column(d):
    cost = 2*d["N"]**2*(d["M"]-d["N"]/3) # for R
    cost += 2*d["N"]**2*(d["M"]-d["N"]/3) # for Q (because it has size m x n)
    return cost

# TODO I think A does not have to be a full rank matrix. In that case, I don't get a full Q.

qr_column = FactorizationKernel(
    Pattern(
        _A,
        lambda d: d["_A"].has_property(properties.FULL_RANK) and \
                  d["_A"].has_property(properties.COLUMN_PANEL) \
                  and not d["_A"].has_property(properties.TRIANGULAR) \
                  and not d["_A"].has_property(properties.DIAGONAL) \
                  and not d["_A"].has_property(properties.ORTHOGONAL) \
                  and not d["_A"].has_property(properties.ORTHOGONAL_COLUMNS)
        ),
    Times(_Q, _R),
    [OutputOperand(_Q, None, ("M", "N"), [properties.FULL_RANK, properties.COLUMN_PANEL, properties.ORTHOGONAL_COLUMNS]),
     OutputOperand(_R, _A, ("N", "N"), [properties.UPPER_TRIANGULAR, properties.SQUARE])
    ],
    cf_qr_column,
    CodeTemplate(),
    CodeTemplate("${type_prefix}geqrf($M, $N, $_A, $ldA, $tau, $work, $lwork, &info); // QR column"),
    CodeTemplate(),
    [SizeArgument("M", _A, "rows"),
     SizeArgument("N", _A, "columns"),
     StrideArgument("ldA", _A, "rows")],
    )

#####################
# Eigendecomposition
#
# syev
# probably no problem

# TODO parts of A are destroyed without being overwritten. Add special thing to specify that.

_A = WildcardSymbol("_A")
_Z = WildcardSymbol("_Z")
_W = WildcardSymbol("_W")
Z = Matrix("Z", (n, n))
W = Matrix("W", (n, n))

def cf_eigen(d):
    cost = 4/3*d["N"]**3 # reduction to tridiagional form
    cost += d["N"]**2 # solve
    cost += 2*d["N"]**3 # backtransformation
    return cost

eigendecomposition = FactorizationKernel(
    Pattern(
        _A,
        lambda d: d["_A"].has_property(properties.SPD) or 
                  (d["_A"].has_property(properties.SYMMETRIC) \
                    and not d["_A"].has_property(properties.DIAGONAL)
                  )
        ),
    Times(_Z, _W, Transpose(_Z)),
    [OutputOperand(_Z, _A, ("N", "N"), [properties.ORTHOGONAL]),
     OutputOperand(_W, None, ("N", "N"), [properties.DIAGONAL, properties.SQUARE])
    ],
    cf_eigen,
    CodeTemplate("""\
                $type* work_$work_id = 0;
                int lwork_$work_id = -1;
                ${type_prefix}syev("V", $uplo, $N, NULL, $ldA, NULL, &work_$work_id, &lwork_$work_id, &info);
                lwork_$work_id = work_$work_id[0];
                work_$work_id = ($type*) malloc(lwork_$work_id);
                """),
    CodeTemplate("""${type_prefix}syev("V", "L", $N, $_A, $ldA, $_W, &work_$work_id, &lwork_$work_id, &info); // eigendecomposition"""),
    CodeTemplate("""free(work_$work_id);\n"""),
    [SizeArgument("N", _A, "rows"),
     StrideArgument("ldA", _A, "rows")]
    )


#####################
# Singular value decomposition
#
# gesvd
# probably no problem

_A = WildcardSymbol("_A")
_U = WildcardSymbol("_U")
_S = WildcardSymbol("_S")
_V = WildcardSymbol("_V")
U = Matrix("U", (m, m))
S = Matrix("S", (m, n))
V = Matrix("V", (n, n))
cf = lambda d: 14*d["M"]*d["N"]**2 + 8*d["N"]**3

# TODO overwriting missing

singular_value = FactorizationKernel(
    Pattern(
        _A,
        # lambda d: d["_A"].has_property(properties.MATRIX)
        # checking those properties is important when applying
        # factorizations to special properties
        lambda d: not d["_A"].has_property(properties.DIAGONAL) and \
                  not d["_A"].has_property(properties.TRIANGULAR) and \
                  not d["_A"].has_property(properties.ORTHOGONAL_COLUMNS) and \
                  not d["_A"].has_property(properties.ORTHOGONAL)
    ),
    Times(_U, _S, _V),
    [OutputOperand(_U, None, ("M", "M"), [properties.ORTHOGONAL]),
     OutputOperand(_S, None, ("M", "N"), [properties.DIAGONAL]),
     OutputOperand(_V, None, ("N", "N"), [properties.ORTHOGONAL])
    ],
    cf,
    CodeTemplate("""\
                $type* work_$work_id = 0;
                int lwork_$work_id = -1;
                ${type_prefix}gesvd("A", "A", $M, $N, NULL, $ldA, NULL, NULL, $ldU, NULL, $ldV, &work_$work_id, &lwork_$work_id, &info);
                lwork_$work_id = work_$work_id[0];
                work_$work_id = ($type*) malloc(lwork_$work_id);
                """),
    CodeTemplate("""${type_prefix}gesvd("A", "A", $M, $N, $_A, $ldA, $_S, $_U, $ldU, $_V, $ldV, &work_$work_id, &lwork_$work_id, &info); // SVD"""),
    CodeTemplate("""free(work_$work_id);\n"""),
    [SizeArgument("M", _A, "rows"),
     SizeArgument("N", _A, "columns"),
     StrideArgument("ldA", _A, "rows"),
     StrideArgument("ldU", _U, "rows"), 
     StrideArgument("ldV", _V, "rows")],
    )

