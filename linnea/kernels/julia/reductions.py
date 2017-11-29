from ..utils.reductions import KernelDescription, \
                               Op1, Op2, \
                               ExpressionKV, PropertyKV, DefaultValueKV, OperatorKV, \
                               KernelType, \
                               OutputOperand

from ..utils.general import SizeArgument, PropertyArgument, StrideArgument, StorageFormatArgument, \
                            InputOperand


from ...code_generation.memory.storage_format import StorageFormat

from ...algebra.properties import Property as properties

from ...algebra.expression import Times, Plus, Transpose, Inverse, Identity, \
                                  ConstantScalar, \
                                  Scalar, Vector, Matrix

n = 10
m = 20
k = 30

###############
# BLAS "0"

# scalar product

alpha = Scalar("alpha")
beta = Scalar("beta")
out = Scalar("out")
cf = lambda d: 1

scalar_product = KernelDescription(
    ExpressionKV(
        None,
        {None: Times(alpha, beta)}
        ),
    [], # variants
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(beta, StorageFormat.full)
    ],
    OutputOperand(out, StorageFormat.full), # return value
    cf, # cost function
    "",
    "$out = $alpha * $beta",
    "",
    [], # Argument objects
    )

# TODO scalar division? -> scalar inversion

# scalar sum

alpha = Scalar("alpha")
beta = Scalar("beta")
out = Scalar("out")
cf = lambda d: 1

scalar_sum = KernelDescription(
    ExpressionKV(
        None,
        {None: Plus(alpha, beta)}
    ),
    [], # variants
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(beta, StorageFormat.full)
    ],
    OutputOperand(out, StorageFormat.full), # return value
    cf, # cost function
    "",
    "$out = $alpha + $beta",
    "",
    [], # Argument objects
    )

###############
# BLAS 1

# DOT

x = Vector("x", (n, 1))
y = Vector("y", (n, 1))
out = Scalar("out")
cf = lambda d: 2*d["N"]

dot = KernelDescription(
    ExpressionKV(
        None,
        {None: Times(Transpose(x),y)}
    ),
    [], # variants
    [InputOperand(x, StorageFormat.full),
     InputOperand(y, StorageFormat.full)
    ],
    OutputOperand(out, StorageFormat.full), # return value
    cf, # cost function
    "",
    "$out = BLAS.dot($N, $x, 1, $y, 1)",
    "",
    [SizeArgument("N", x, "rows")], # Argument objects
    )

# SCAL

x = Vector("x", (n, 1))
alpha = Scalar("alpha")
cf = lambda d: d["N"]

scal = KernelDescription(
    ExpressionKV(
        None,
        {None: Times(alpha, x)}
    ),
    [], # variants
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(x, StorageFormat.full)
    ],
    OutputOperand(x, StorageFormat.full), # return value
    cf, # cost function
    "",
    "scal!($N, $alpha, $x, 1)",
    "",
    [SizeArgument("N", x, "rows")], # Argument objects
    )

# AXPY

x = Vector("x", (n, 1))
y = Vector("y", (n, 1))
alpha = Scalar("alpha")
cf = lambda d: d["N"]

axpy = KernelDescription(
    ExpressionKV(
        None,
        {None: Plus(Times(alpha, x), y)}),
        # {None: Plus(Times(alpha, Op1(x)), Op2(y))}),
    [
        DefaultValueKV(alpha, [ConstantScalar(1.0)]),
        # OperatorKV("", {"N": Identity, "T": Transpose}, Op1),
        # OperatorKV("", {"N": Identity, "T": Transpose}, Op2),
    ],
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(x, StorageFormat.full),
     InputOperand(y, StorageFormat.full)
    ],
    OutputOperand(y, StorageFormat.full), # return value
    cf, # cost function
    "",
    "axpy!($alpha, $x, $y) # vectors",
    "",
    [SizeArgument("N", x, "rows")], # Argument objects
    )

###############
# BLAS 2

# GER (outer product)

A = Matrix("A", (m, n))
x = Vector("x", (m, 1))
y = Vector("y", (n, 1))
alpha = Scalar("alpha")
cf = lambda d: 2*d["N"]*d["M"]

ger = KernelDescription(
    ExpressionKV(
        None,
        {None: Plus(Times(alpha, x, Transpose(y)), A)}
    ),
    [
        DefaultValueKV(alpha, [ConstantScalar(1.0)])
    ],
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(x, StorageFormat.full),
     InputOperand(y, StorageFormat.full),
     InputOperand(A, StorageFormat.full)
    ],
    OutputOperand(A, StorageFormat.full), # return value
    cf, # cost function
    "",
    "ger!($alpha, $x, $y, $A)",
    "",
    [SizeArgument("M", x, "rows"),
     SizeArgument("N", y, "rows")], # Argument objects
    )

A = Matrix("A", (m, n))
x = Vector("x", (m, 1))
y = Vector("y", (n, 1))
alpha = Scalar("alpha")
cf = lambda d: 2*d["N"]*d["M"]

ger_alt = KernelDescription(
    ExpressionKV(
        None,
        {None: Times(alpha, x, Transpose(y))}
    ),
    [
        DefaultValueKV(alpha, [ConstantScalar(1.0)])
    ],
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(x, StorageFormat.full),
     InputOperand(y, StorageFormat.full)
    ],
    OutputOperand(A, StorageFormat.full), # return value
    cf, # cost function
    "",
    # If A is not allocated, use "$A = zeros($type, ($M, $N))\n"
    "fill!($A, 0.0)\nger!($alpha, $x, $y, $A)",
    "",
    [SizeArgument("M", x, "rows"),
     SizeArgument("N", y, "rows")], # Argument objects
    )


# SYR (outer product)

A = Matrix("A", (m, n))
x = Vector("x", (m, 1))
alpha = Scalar("alpha")
cf = lambda d: d["M"]**2

syr = KernelDescription(
    ExpressionKV(
        None,
        {None: Plus(Times(alpha, x, Transpose(x)), A)}
    ),
    [
        DefaultValueKV(alpha, [ConstantScalar(1.0)])
    ],
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(x, StorageFormat.full),
     InputOperand(A, StorageFormat.full)
    ],
    OutputOperand(A, StorageFormat.symmetric_lower_triangular), # return value
    cf, # cost function
    "",
    "syr!('L', $alpha, $x, $A)",
    "",
    [SizeArgument("M", x, "rows")], # Argument objects
    )

A = Matrix("A", (m, n))
x = Vector("x", (m, 1))
alpha = Scalar("alpha")
cf = lambda d: d["M"]**2

syr_alt = KernelDescription(
    ExpressionKV(
        None,
        {None: Times(alpha, x, Transpose(x))}
    ),
    [
        DefaultValueKV(alpha, [ConstantScalar(1.0)])
    ],
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(x, StorageFormat.full),
    ],
    OutputOperand(A, StorageFormat.symmetric_lower_triangular), # return value
    cf, # cost function
    "",
    # If A is not allocated, use "$A = zeros($type, ($M, $N))\n"
    "fill!($A, 0.0)\nsyr!('L', $alpha, $x, $A)",
    "",
    [SizeArgument("M", x, "rows")], # Argument objects
    )


# GEMV

A = Matrix("A", (m, n))
x = Vector("x", (n, 1))
y = Vector("y", (m, 1))
alpha = Scalar("alpha")
beta = Scalar("beta")
cf = lambda d: 2*d["N"]*d["M"]

"""
TODO also use the following?
- gemv(tA, alpha, A, x)
- gemv(tA, A, x)
Compare performance.
"""

gemv = KernelDescription(
    ExpressionKV(
        None,
        {None: Plus(Times(alpha, Op1(A), x), Times(beta, y))}
    ),
    [
        OperatorKV("transA", {"N": Identity, "T": Transpose}, Op1),
        DefaultValueKV(alpha, [ConstantScalar(1.0)]),
        DefaultValueKV(beta, [ConstantScalar(0.0), ConstantScalar(1.0)])
    ],
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(A, StorageFormat.full),
     InputOperand(x, StorageFormat.full),
     InputOperand(beta, StorageFormat.full),
     InputOperand(y, StorageFormat.full),
    ],
    OutputOperand(y, StorageFormat.full), # return value
    cf, # cost function
    "",
    "gemv!($transA, $alpha, $A, $x, $beta, $y)",
    "",
    [SizeArgument("M", A, "rows"),
     SizeArgument("N", A, "columns")], # Argument objects
    [KernelType.identity, KernelType.transpose]
    )


# SYMV

A = Matrix("A", (m, n))
A.set_property(properties.SYMMETRIC)
x = Vector("x", (n, 1))
y = Vector("y", (m, 1))
alpha = Scalar("alpha")
beta = Scalar("beta")
cf = lambda d: d["M"]**2

symv = KernelDescription(
    ExpressionKV(
        None,
        {None: Plus(Times(alpha, A, x), Times(beta, y))}
    ),
    [
        DefaultValueKV(alpha, [ConstantScalar(1.0)]),
        DefaultValueKV(beta, [ConstantScalar(0.0), ConstantScalar(1.0)])
    ],
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(A, StorageFormat.symmetric_triangular),
     InputOperand(x, StorageFormat.full),
     InputOperand(beta, StorageFormat.full),
     InputOperand(y, StorageFormat.full),
    ],
    OutputOperand(y, StorageFormat.full), # return value
    cf, # cost function
    "",
    "symv!($uplo, $alpha, $A, $x, $beta, $y)",
    "",
    [SizeArgument("M", A, "rows"),
     StorageFormatArgument("uplo", A, StorageFormat.symmetric_lower_triangular, ["L", "U"]),], # Argument objects
    [KernelType.identity, KernelType.transpose]
    )


# TRMV

"""
TODO: Julia documentation wrong
"""

A = Matrix("A", (n, n))
A.set_property(properties.SQUARE)
x = Vector("x", (n, 1))
cf = lambda d: d["N"]**2

trmv = KernelDescription(
    ExpressionKV(
        None,
        {None: Times(Op1(A), x)}
    ),
    [
        OperatorKV("transA", {"N": Identity, "T": Transpose}, Op1),
        PropertyKV("uplo", {"U": properties.UPPER_TRIANGULAR, "L": properties.LOWER_TRIANGULAR}, A)
    ],
    [InputOperand(A, StorageFormat.triangular_udiag_opt),
     InputOperand(x, StorageFormat.full),
    ],
    OutputOperand(x, StorageFormat.full), # return value
    cf, # cost function
    "",
    "trmv!($uplo, $transA, $diag, $A, $x)",
    "",
    [SizeArgument("N", A, "rows"),
     PropertyArgument("diag", A, properties.UNIT_DIAGONAL, ["U", "N"])], # Argument objects
    [KernelType.identity, KernelType.transpose]
    )


# TRSV

A = Matrix("A", (n, n))
A.set_property(properties.SQUARE)
x = Vector("x", (n, 1))
cf = lambda d: d["N"]**2

"""
TODO also use the following?
- trsv(ul, tA, dA, A, b)
Compare performance.
"""

trsv = KernelDescription(
    ExpressionKV(
        None,
        {None: Times(Op1(Inverse(A)), x)}
    ),
    [
        OperatorKV("transA", {"N": Identity, "T": Transpose}, Op1),
        PropertyKV("uplo", {"U": properties.UPPER_TRIANGULAR, "L": properties.LOWER_TRIANGULAR}, A)
    ],
    [InputOperand(A, StorageFormat.triangular_udiag_opt),
     InputOperand(x, StorageFormat.full),
    ],
    OutputOperand(x, StorageFormat.full), # return value
    cf, # cost function
    "",
    "trsv!($uplo, $transA, $diag, $A, $x)",
    "",
    [SizeArgument("N", A, "rows"),
     PropertyArgument("diag", A, properties.UNIT_DIAGONAL, ["U", "N"])], # Argument objects
    [KernelType.identity, KernelType.transpose]
    )


###############
# BLAS 3

# GEMM

A = Matrix("A", (m, k))
B = Matrix("B", (k, n))
C = Matrix("C", (m, n))
alpha = Scalar("alpha")
beta = Scalar("beta")
cf = lambda d: 2*d["N"]*d["M"]*d["K"]

"""
TODO also use the following?
- gemm(tA, tB, alpha, A, B)
- gemm(tA, tB, A, B)
Compare performance.
"""

gemm = KernelDescription(
    ExpressionKV(
        None,
        {None: Plus(Times(alpha, Op1(A), Op2(B)), Times(beta, C))}
    ),
    [
        OperatorKV("transA", {"N": Identity, "T": Transpose}, Op1),
        OperatorKV("transB", {"N": Identity, "T": Transpose}, Op2),
        DefaultValueKV(alpha, [ConstantScalar(1.0)]),
        DefaultValueKV(beta, [ConstantScalar(0.0), ConstantScalar(1.0)])
    ],
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(A, StorageFormat.full),
     InputOperand(B, StorageFormat.full),
     InputOperand(beta, StorageFormat.full),
     InputOperand(C, StorageFormat.full),
    ],
    OutputOperand(C, StorageFormat.full), # return value
    cf, # cost function
    "",
    "gemm!($transA, $transB, $alpha, $A, $B, $beta, $C)",
    "",
    [SizeArgument("M", Op1(A), "rows"),
     SizeArgument("N", Op2(B), "columns"),
     SizeArgument("K", Op1(A), "columns")], # Argument objects
    )


# SYMM

A = Matrix("A", (m, k))
A.set_property(properties.SYMMETRIC)
B = Matrix("B", (k, n))
C = Matrix("C", (m, n))
alpha = Scalar("alpha")
beta = Scalar("beta")
cf = lambda d: d["M"]*d["N"]*d["K"]
"""
The names for the sizes in this cost function do not make a lot of sense.
However, it is not possible to do something more consistent because of the
side argument. A more intuitive way to describe the cost would be n^2*m, where
n is the size of matrix A, and m is the "other" size.
"""

symm = KernelDescription(
    ExpressionKV(
        "side",
        {"L": Plus(Times(alpha, A, B), Times(beta, C)),
         "R": Plus(Times(alpha, B, A), Times(beta, C))}
    ),
    [
        DefaultValueKV(alpha, [ConstantScalar(1.0)]),
        DefaultValueKV(beta, [ConstantScalar(0.0), ConstantScalar(1.0)])
    ],
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(A, StorageFormat.symmetric_triangular),
     InputOperand(B, StorageFormat.full),
     InputOperand(beta, StorageFormat.full),
     InputOperand(C, StorageFormat.full),
    ],
    OutputOperand(C, StorageFormat.full), # return value
    cf, # cost function
    "",
    "symm!($side, $uplo, $alpha, $A, $B, $beta, $C)",
    "",
    [SizeArgument("M", B, "rows"),
     SizeArgument("N", B, "columns"),
     SizeArgument("K", A, "rows"),
     StorageFormatArgument("uplo", A, StorageFormat.symmetric_lower_triangular, ["L", "U"])], # Argument objects
    )



# SYRK

A = Matrix("A", (m, k))
C = Matrix("C", (m, n))
C.set_property(properties.SYMMETRIC)
alpha = Scalar("alpha")
beta = Scalar("beta")
cf = lambda d: d["K"]*(d["M"]**2)


syrk = KernelDescription(
    ExpressionKV(
        None,
        {None: Plus(Times(alpha, Op1(A), Op1(Transpose(A))), Times(beta, C))}
    ),
    [
        OperatorKV("transA", {"N": Identity, "T": Transpose}, Op1),
        DefaultValueKV(alpha, [ConstantScalar(1.0)]),
        DefaultValueKV(beta, [ConstantScalar(0.0), ConstantScalar(1.0)])
    ],
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(A, StorageFormat.full),
     InputOperand(beta, StorageFormat.full),
     InputOperand(C, StorageFormat.symmetric_triangular),
    ],
    OutputOperand(C, StorageFormat.symmetric_triangular_out), # return value
    cf, # cost function
    "",
    "syrk!($uplo, $transA, $alpha, $A, $beta, $C)",
    "",
    [SizeArgument("M", Op1(A), "rows"),
     SizeArgument("K", Op1(A), "columns"),
     StorageFormatArgument("uplo", C, StorageFormat.symmetric_lower_triangular, ["L", "U"])], # Argument objects
    )


# TRMM

A = Matrix("A", (m, m))
A.set_property(properties.SQUARE)
B = Matrix("B", (m, n))
alpha = Scalar("alpha")
cf = lambda d: d["M"]*d["N"]*d["K"]
"""
The names for the sizes in this cost function do not make a lot of sense.
However, it is not possible to do something more consistent because of the
side argument. A more intuitive way to describe the cost would be n^2*m, where
n is the size of matrix A, and m is the "other" size.
"""


trmm = KernelDescription(
    ExpressionKV(
        "side",
        {"L": Times(alpha, Op1(A), B),
         "R": Times(alpha, B, Op1(A))}
    ),
    [
        OperatorKV("transA", {"N": Identity, "T": Transpose}, Op1),
        DefaultValueKV(alpha, [ConstantScalar(1.0)]),
        PropertyKV("uplo", {"U": properties.UPPER_TRIANGULAR, "L": properties.LOWER_TRIANGULAR}, A)
    ],
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(A, StorageFormat.triangular_udiag_opt),
     InputOperand(B, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    "",
    "trmm!($side, $uplo, $transA, $diag, $alpha, $A, $B)",
    "",
    [SizeArgument("M", B, "rows"),
     SizeArgument("N", B, "columns"),
     SizeArgument("K", A, "rows"),
     PropertyArgument("diag", A, properties.UNIT_DIAGONAL, ["U", "N"])], # Argument objects
    [KernelType.identity, KernelType.transpose]
    )


# TRSM

A = Matrix("A", (m, m))
A.set_property(properties.SQUARE)
B = Matrix("B", (m, n))
alpha = Scalar("alpha")
cf = lambda d: d["M"]*d["N"]*d["K"]
"""
The names for the sizes in this cost function do not make a lot of sense.
However, it is not possible to do something more consistent because of the
side argument. A more intuitive way to describe the cost would be n^2*m, where
n is the size of matrix A, and m is the "other" size.
"""

"""
TODO also use the following?
- trsm(side, ul, tA, dA, alpha, A, B)
Compare performance.
"""

trsm = KernelDescription(
    ExpressionKV(
        "side",
        {"L": Times(alpha, Op1(Inverse(A)), B),
         "R": Times(alpha, B, Op1(Inverse(A)))}
    ),
    [
        OperatorKV("transA", {"N": Identity, "T": Transpose}, Op1),
        DefaultValueKV(alpha, [ConstantScalar(1.0)]),
        PropertyKV("uplo", {"U": properties.UPPER_TRIANGULAR, "L": properties.LOWER_TRIANGULAR}, A)
    ],
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(A, StorageFormat.triangular_udiag_opt),
     InputOperand(B, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    "",
    "trsm!($side, $uplo, $transA, $diag, $alpha, $A, $B)",
    "",
    [SizeArgument("M", B, "rows"),
     SizeArgument("N", B, "columns"),
     SizeArgument("K", A, "rows"),
     PropertyArgument("diag", A, properties.UNIT_DIAGONAL, ["U", "N"])], # Argument objects
    [KernelType.identity, KernelType.transpose]
    )

###############
# LAPACK

# SCAL (scalar * matrix) (equivalent to LASCL in LAPACK)

X = Matrix("X", (m, n))
alpha = Scalar("alpha")
cf = lambda d: d["M"]*d["N"]

lascl = KernelDescription(
    ExpressionKV(
        None,
        {None: Times(alpha, X)}
    ),
    [], # variants
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(X, StorageFormat.full),
    ],
    OutputOperand(X, StorageFormat.full), # return value
    cf, # cost function
    "",
    "scal!($MN, $alpha, $X, 1)",
    "",
    [SizeArgument("M", X, "rows"),
     SizeArgument("N", X, "columns"),
     SizeArgument("MN", X, "entries")], # Argument objects
    )


# getri (matrix inversion)

A = Matrix("A", (n, n))
B = Matrix("B", (n, n))
A.set_property(properties.SQUARE)
cf = lambda d: 2*d["N"]**3

getri = KernelDescription(
    ExpressionKV(
        None,
        {None: Inverse(A)}
    ),
    [], # variants
    [InputOperand(A, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    "",
    "$B = inv($A)",
    "",
    [SizeArgument("N", A, "columns")], # Argument objects
    [KernelType.identity, KernelType.transpose]
    )

# trtri (triangular matrix inversion)

A = Matrix("A", (n, n))
B = Matrix("B", (n, n))
A.set_property(properties.SQUARE)
A.set_property(properties.TRIANGULAR)
cf = lambda d: (d["N"]**3)/3

# TODO this can not be used at the moment because of StorageFormat.as_overwritten
# and storage format conversions

trtri = KernelDescription(
    ExpressionKV(
        None,
        {None: Inverse(A)}
    ),
    [
        PropertyKV("uplo", {"U": properties.UPPER_TRIANGULAR, "L": properties.LOWER_TRIANGULAR}, A)
    ], # variants
    [InputOperand(A, StorageFormat.triangular_udiag_opt),
    ],
    OutputOperand(A, StorageFormat.as_overwritten), # return value
    cf, # cost function
    "",
    "Base.LinAlg.LAPACK.trtri!($uplo, $diag, $A)",
    "",
    [SizeArgument("N", A, "rows"),
     PropertyArgument("diag", A, properties.UNIT_DIAGONAL, ["U", "N"])], # Argument objects
    [KernelType.identity, KernelType.transpose]
    )

# trtri = KernelDescription(
#     ExpressionKV(
#         None,
#         {None: Inverse(A)}
#     ),
#     [], # variants
#     [InputOperand(A, StorageFormat.full),
#     ],
#     OutputOperand(B, StorageFormat.full), # return value
#     cf, # cost function
#     "",
#     "$B = inv($A)",
#     "",
#     [SizeArgument("N", A, "rows")], # Argument objects
#     [KernelType.identity, KernelType.transpose]
#     )


# POSV

A = Matrix("A", (m, m))
A.set_property(properties.SQUARE)
A.set_property(properties.SPD)
B = Matrix("B", (m, n))
cf = lambda d: (d["N"]**3)/3 + 2*(d["N"]**2)*d["M"]

"""
TODO problem: both A and B are overwritten, but it's not possible to express that here
alternative
Base.LinAlg.LAPACK.potrf!('L', A)
Base.LinAlg.LAPACK.potrs!('L', A, B)
"""

posv = KernelDescription(
    ExpressionKV(
        None,
        {None: Times(Inverse(A), B)}
    ),
    [],  # variants
    [InputOperand(A, StorageFormat.symmetric_triangular),
     InputOperand(B, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    "",
    "Base.LinAlg.LAPACK.posv!('L', $A, $B)",
    "",
    [SizeArgument("M", B, "rows"),
     SizeArgument("N", B, "columns")], # Argument objects
    [KernelType.identity, KernelType.transpose]
    )


# POSV right

A = Matrix("A", (m, m))
A.set_property(properties.SQUARE)
A.set_property(properties.SPD)
B = Matrix("B", (n, m))
cf = lambda d: (d["N"]**3)/3 + 2*(d["N"]**2)*d["M"]

"""
TODO problem: both A and B are overwritten, but it's not possible to express that here
alternative
Base.LinAlg.LAPACK.potrf!('L', A)
Base.LinAlg.LAPACK.potrs!('L', A, B)
"""

posvr = KernelDescription(
    ExpressionKV(
        None,
        {None: Times(B, Inverse(A))}
    ),
    [],  # variants
    [InputOperand(A, StorageFormat.symmetric_triangular),
     InputOperand(B, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    "",
    "$B = $B'\nBase.LinAlg.LAPACK.posv!('L', $A, $B)\n$B = $B'",
    "",
    [SizeArgument("M", B, "columns"),
     SizeArgument("N", B, "rows")], # Argument objects
    [KernelType.identity, KernelType.transpose]
    )


# sysv

A = Matrix("A", (m, m))
A.set_property(properties.SQUARE)
A.set_property(properties.SYMMETRIC)
B = Matrix("B", (m, n))
cf = lambda d: (d["N"]**3)/3 + 2*(d["N"]**2)*d["M"]

"""
TODO problem: both A and B are overwritten, but it's not possible to express that here
alternative
(A, ipiv) = Base.LinAlg.LAPACK.sytrf!('L', A)
Base.LinAlg.LAPACK.sytrs!('L', A, ipiv, B)
TODO For whatever reason, sytrs is very slow. Investigate.
"""

sysv = KernelDescription(
    ExpressionKV(
        None,
        {None: Times(Inverse(A), B)}
    ),
    [],  # variants
    [InputOperand(A, StorageFormat.symmetric_triangular),
     InputOperand(B, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    "",
    "Base.LinAlg.LAPACK.sysv!('L', $A, $B)",
    "",
    [SizeArgument("M", B, "rows"),
     SizeArgument("N", B, "columns")], # Argument objects
    [KernelType.identity, KernelType.transpose]
    )


# sysv right

A = Matrix("A", (m, m))
A.set_property(properties.SQUARE)
A.set_property(properties.SYMMETRIC)
B = Matrix("B", (n, m))
cf = lambda d: (d["N"]**3)/3 + 2*(d["N"]**2)*d["M"]

"""
TODO problem: both A and B are overwritten, but it's not possible to express that here
alternative
(A, ipiv) = Base.LinAlg.LAPACK.sytrf!('L', A)
Base.LinAlg.LAPACK.sytrs!('L', A, ipiv, B)
TODO For whatever reason, sytrs is very slow. Investigate.
"""

sysvr = KernelDescription(
    ExpressionKV(
        None,
        {None: Times(B, Inverse(A))}
    ),
    [],  # variants
    [InputOperand(A, StorageFormat.symmetric_triangular),
     InputOperand(B, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    "",
    "$B = $B'\nBase.LinAlg.LAPACK.sysv!('L', $A, $B)\n$B = $B'",
    "",
    [SizeArgument("M", B, "columns"),
     SizeArgument("N", B, "rows")], # Argument objects
    [KernelType.identity, KernelType.transpose]
    )

# gesv

A = Matrix("A", (m, m))
A.set_property(properties.SQUARE)
B = Matrix("B", (m, n))
cf = lambda d: 2*(d["N"]**3)/3 + 2*(d["N"]**2)*d["M"]

"""
TODO problem: both A and B are overwritten, but it's not possible to express that here
alternative
Base.LinAlg.LAPACK.gesv!($A, $B)
"""

gesv = KernelDescription(
    ExpressionKV(
        None,
        {None: Times(Op1(Inverse(A)), B)}
    ),
    [
        OperatorKV("transA", {"N": Identity, "T": Transpose}, Op1),
    ],  # variants
    [InputOperand(A, StorageFormat.full),
     InputOperand(B, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    "",
    "($A, ipiv, info) = Base.LinAlg.LAPACK.getrf!($A)\nBase.LinAlg.LAPACK.getrs!($transA, $A, ipiv, $B)",
    "",
    [SizeArgument("M", B, "rows"),
     SizeArgument("N", B, "columns")], # Argument objects
    [KernelType.identity, KernelType.transpose]
    )


# gesv right

A = Matrix("A", (m, m))
A.set_property(properties.SQUARE)
B = Matrix("B", (n, m))
cf = lambda d: 2*(d["N"]**3)/3 + 2*(d["N"]**2)*d["M"]

"""
TODO problem: both A and B are overwritten, but it's not possible to express that here
alternative
Base.LinAlg.LAPACK.gesv!($A, $B)
"""

gesvr = KernelDescription(
    ExpressionKV(
        None,
        {None: Times(B, Inverse(A))}
    ),
    [],  # variants
    [InputOperand(A, StorageFormat.full),
     InputOperand(B, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    "",
    "$B = $B/lufact!($A)",
    "",
    [SizeArgument("M", B, "columns"),
     SizeArgument("N", B, "rows")], # Argument objects
    [KernelType.identity, KernelType.transpose]
    )

# gesv right transpose

A = Matrix("A", (m, m))
A.set_property(properties.SQUARE)
B = Matrix("B", (n, m))
cf = lambda d: 2*(d["N"]**3)/3 + 2*(d["N"]**2)*d["M"]

"""
TODO problem: both A and B are overwritten, but it's not possible to express that here
alternative
Base.LinAlg.LAPACK.gesv!($A, $B)
"""

gesvrt = KernelDescription(
    ExpressionKV(
        None,
        {None: Times(B, Inverse(A))}
    ),
    [],  # variants
    [InputOperand(A, StorageFormat.full),
     InputOperand(B, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    "",
    "$B = $B/lufact!($A')",
    "",
    [SizeArgument("M", B, "columns"),
     SizeArgument("N", B, "rows")], # Argument objects
    [KernelType.identity, KernelType.transpose]
    )

# diaginv (diagonal matrix inversion)

A = Matrix("A", (n, n))
A.set_property(properties.SQUARE)
A.set_property(properties.DIAGONAL)
cf = lambda d: d["N"]

diaginv = KernelDescription(
    ExpressionKV(
        None,
        {None: Inverse(A)}
    ),
    [], # variants
    [InputOperand(A, StorageFormat.diagonal_vector),
    ],
    OutputOperand(A, StorageFormat.diagonal_vector), # return value
    cf, # cost function
    "",
    "$A = 1./$A",
    "",
    [SizeArgument("N", A, "rows")], # Argument objects
    [KernelType.identity, KernelType.transpose]
    )

# matrix transposition

A = Matrix("A", (m, n))
B = Matrix("B", (n, m))
cf = lambda d: 0

transpose = KernelDescription(
    ExpressionKV(
        None,
        {None: Transpose(A)}
    ),
    [], # variants
    [InputOperand(A, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    "",
    # $B = Array{$type}($N, $M)
    """\
    transpose!($B, $A)\
    """,
    "",
    [SizeArgument("M", A, "rows"),
     SizeArgument("N", A, "columns")], # Argument objects
    # [KernelType.identity, KernelType.transpose]
    )

# vector transposition (totally fake) (this is kind of tricky. Since it doesn't actually do anything, overwriting is never a problem. Solution: Just use different symbol for output?)

x = Vector("x", (n, 1))
cf = lambda d: 0

transpose_vector = KernelDescription(
    ExpressionKV(
        None,
        {None: Transpose(x)}
    ),
    [], # variants
    [InputOperand(x, StorageFormat.full),
    ],
    OutputOperand(x, StorageFormat.full), # return value
    cf, # cost function
    "",
    "# Transposing vector $x (no operation);",
    "",
    [], # Argument objects
    # [KernelType.identity, KernelType.transpose]
    )

# matrix sum

alpha = Scalar("alpha")
A = Matrix("A", (m, n))
B = Matrix("B", (m, n))
cf = lambda d: 2*d["N"]*d["M"]

matrix_sum = KernelDescription(
    ExpressionKV(
        None,
        {None: Plus(Times(alpha, A), B)} # a*x+y
    ),
    [
        DefaultValueKV(alpha, [ConstantScalar(1.0)])
    ],
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(A, StorageFormat.full),
     InputOperand(B, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    "",
    "axpy!($alpha, $A, $B) # matrices",
    "",
    [SizeArgument("N", A, "columns"),
     SizeArgument("M", A, "rows")], # Argument objects
    )

# scalar * diagonal

X = Matrix("X", (m, n))
X.set_property(properties.DIAGONAL)
alpha = Scalar("alpha")
cf = lambda d: min(d["M"], d["N"])

diagscal = KernelDescription(
    ExpressionKV(
        None,
        {None: Times(alpha, X)}
    ),
    [], # variants
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(X, StorageFormat.diagonal_vector),
    ],
    OutputOperand(X, StorageFormat.diagonal_vector), # return value
    cf, # cost function
    "",
    "scal!(min($M, $N), $alpha, $X, 1)",
    "",
    [SizeArgument("M", X, "rows"),
     SizeArgument("N", X, "columns")], # Argument objects
    )

# diagonal * diagonal

"""
TODO: For "outer products", parts of the diagonal have to be
filled with zeros. For now, this is done with a fairly ugly hack. Can it be done simpler?

TODO: what about transA, transB? Since diagonal matrices don't have to be
square, it could happen. I think it wouldn't even change the code.

Remark: Overwriting is not possible because both A and B could have a different
length than C (because diagonal matrices don't have to be square). To make
overwriting possible in some cases, one could add a variant where one operand
has to be square. That's perhaps a fairly common case. With a reasonable method
for the selection of kernels, this kernel will be preferred because it is more
specific, even if the number of FLOPS is the same.
"""

A = Matrix("A", (m, k))
A.set_property(properties.DIAGONAL)
B = Matrix("B", (k, n))
B.set_property(properties.DIAGONAL)
C = Matrix("C", (m, n))
cf = lambda d: min(d["M"], d["K"], d["N"])

diagdiagmul = KernelDescription(
    ExpressionKV(
        None,
        {None: Times(Op1(A), Op2(B))}
    ),
    [
        OperatorKV("", {"N": Identity, "T": Transpose}, Op1),
        OperatorKV("", {"N": Identity, "T": Transpose}, Op2),
    ], # variants
    [InputOperand(A, StorageFormat.diagonal_vector),
     InputOperand(B, StorageFormat.diagonal_vector),
    ],
    OutputOperand(C, StorageFormat.diagonal_vector), # return value
    cf, # cost function
    "",
    "for i = 1:min(length($A), length($B)); $C[i] = $A[i] * $B[i]; end; if max(length($A), length($B)) < length($C) $C[min(length($A), length($B))+1:end] = 0; end;",
    "",
    [SizeArgument("M", A, "rows"),
     SizeArgument("K", A, "columns"),
     SizeArgument("N", B, "columns")], # Argument objects
    )


# diagonal * inv(diagonal)
# inv(diagonal) * diagonal

"""
TODO: what about transB? Since B doesn't have to be
square, it could happen.
"""

A = Matrix("A", (m, m))
A.set_property(properties.DIAGONAL)
A.set_property(properties.SQUARE)
B = Matrix("B", (m, n))
B.set_property(properties.DIAGONAL)
cf = lambda d: min(d["M"], d["N"])

diagdiagsolve = KernelDescription(
    ExpressionKV(
        "side",
        {"L": Times(Inverse(A), B),
         "R": Times(B, Inverse(A))}
    ),
    [], # variants
    [InputOperand(A, StorageFormat.diagonal_vector),
     InputOperand(B, StorageFormat.diagonal_vector),
    ],
    OutputOperand(B, StorageFormat.diagonal_vector), # return value
    cf, # cost function
    "",
    "for i = 1:length($B); $B[i] /= $A[i]; end;",
    "",
    [SizeArgument("M", A, "rows"),
     SizeArgument("N", B, "columns")], # Argument objects
    [KernelType.identity, KernelType.transpose]
    )


# full * inv(diag)

# scaling columns

"""
TODO: what about transB?
"""

A = Matrix("A", (m, m))
A.set_property(properties.DIAGONAL)
A.set_property(properties.SQUARE)
B = Matrix("B", (m, n))
cf = lambda d: d["M"]*d["N"]

diagsmr = KernelDescription(
    ExpressionKV(
        None,
        {None: Times(B, Inverse(A))}
    ),
    [], # variants
    [InputOperand(A, StorageFormat.diagonal_vector),
     InputOperand(B, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    "",
    "x = 1./$A; for i = 1:size($B, 2); for j=1:size($B, 1); $B[j,i] *= x[i]; end; end;", # faster with current Julia version
    # "for i = 1:size($B, 2); $B[:,i] /= $A[i]; end;",
    "",
    [SizeArgument("M", A, "rows"),
     SizeArgument("N", B, "columns")], # Argument objects
    [KernelType.identity, KernelType.transpose]
    )


# inv(diag) * full

# scaling rows

"""
TODO: what about transB?
"""

A = Matrix("A", (m, m))
A.set_property(properties.DIAGONAL)
A.set_property(properties.SQUARE)
B = Matrix("B", (m, n))
cf = lambda d: d["M"]*d["N"]

diagsml = KernelDescription(
    ExpressionKV(
        None,
        {None: Times(Inverse(A), B)}
    ),
    [], # variants
    [InputOperand(A, StorageFormat.diagonal_vector),
     InputOperand(B, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    "",
    "x = 1./$A; for i = 1:size($B, 1); for j=1:size($B, 2); $B[i,j] *= x[i]; end; end;", # faster with current Julia version
    # "for i = 1:size($B, 1); $B[i,:] /= $A[i]; end;",
    "",
    [SizeArgument("M", A, "rows"),
     SizeArgument("N", B, "columns")], # Argument objects
    [KernelType.identity, KernelType.transpose]
    )


# inv(diag) * vector

"""
TODO: does the transpose variant get simplified correctly? (i.e. A is symmetric)
"""

A = Matrix("A", (m, m))
A.set_property(properties.DIAGONAL)
A.set_property(properties.SQUARE)
x = Vector("x", (m, 1))
cf = lambda d: d["M"]

diagsv = KernelDescription(
    ExpressionKV(
        None,
        {None: Times(Inverse(A), x)}
    ),
    [], # variants
    [InputOperand(A, StorageFormat.diagonal_vector),
     InputOperand(x, StorageFormat.full),
    ],
    OutputOperand(x, StorageFormat.full), # return value
    cf, # cost function
    "",
    "$x ./= $A",
    "",
    [SizeArgument("M", A, "rows")], # Argument objects
    [KernelType.identity, KernelType.transpose]
    )


# full * diag (square)

# scaling columns

"""
TODO: what about transB?

scale!(A, b) is deprecated
"""

A = Matrix("A", (m, m))
A.set_property(properties.DIAGONAL)
A.set_property(properties.SQUARE)
B = Matrix("B", (m, n))
cf = lambda d: d["M"]*d["N"]

diagmmr = KernelDescription(
    ExpressionKV(
        None,
        {None: Times(B, A)}
    ),
    [], # variants
    [InputOperand(A, StorageFormat.diagonal_vector),
     InputOperand(B, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    "",
    "for i = 1:size($B, 2); for j=1:size($B, 1); $B[j,i] *= $A[i]; end; end;", # faster with current Julia version
    # "for i = 1:size($B, 2); $B[:,i] *= $A[i]; end;",
    "",
    [SizeArgument("M", A, "rows"),
     SizeArgument("N", B, "columns")], # Argument objects
    [KernelType.identity, KernelType.transpose]
    )


# diag (square) * full

# scaling rows

"""
TODO: what about transB?

scale!(b, A) is deprecated
"""

A = Matrix("A", (m, m))
A.set_property(properties.DIAGONAL)
A.set_property(properties.SQUARE)
B = Matrix("B", (m, n))
cf = lambda d: d["M"]*d["N"]

diagmml = KernelDescription(
    ExpressionKV(
        None,
        {None: Times(A, B)}
    ),
    [], # variants
    [InputOperand(A, StorageFormat.diagonal_vector),
     InputOperand(B, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    "",
    "for i = 1:size($B, 1); for j=1:size($B, 2); $B[i,j] *= $A[i]; end; end;", # faster with current Julia version
    # "for i = 1:size($B, 1); $B[i,:] *= $A[i]; end;",
    "",
    [SizeArgument("M", A, "rows"),
     SizeArgument("N", B, "columns")], # Argument objects
    [KernelType.identity, KernelType.transpose]
    )


# diag (square) * vector

"""
TODO: does the transpose variant get simplified correctly? (i.e. A is symmetric)
"""

A = Matrix("A", (m, m))
A.set_property(properties.DIAGONAL)
A.set_property(properties.SQUARE)
x = Vector("x", (m, 1))
cf = lambda d: d["M"]

diagmv = KernelDescription(
    ExpressionKV(
        None,
        {None: Times(A, x)}
    ),
    [], # variants
    [InputOperand(A, StorageFormat.diagonal_vector),
     InputOperand(x, StorageFormat.full),
    ],
    OutputOperand(x, StorageFormat.full), # return value
    cf, # cost function
    "",
    "$x .*= $A",
    "",
    [SizeArgument("M", A, "rows")], # Argument objects
    [KernelType.identity, KernelType.transpose]
    )


# matrix and vector
# A[:,F[:p]] # A*P
# A[:,invperm(p)] # A*P^T
# A[invperm(p),:] # P*A
# A[p,:] # P^T*A
# perumation*permutation
# what about diagonal A? could be almost the same as vector

