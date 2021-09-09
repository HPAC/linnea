from ..utils.reductions import KernelDescription, \
                               Op1, Op2, \
                               Operation, OperationKV, \
                               PropertyKV, DefaultValueKV, OperatorKV, \
                               KernelOption, \
                               OutputOperand

from ..utils.general import SizeArgument, PropertyArgument, StrideArgument, \
                            StorageFormatArgument, InputOperand


from ...code_generation.memory.storage_format import StorageFormat

from ...algebra.properties import Property

from ...algebra.expression import Times, Plus, \
                                  Transpose, Inverse, InverseTranspose, \
                                  Identity, \
                                  ConstantScalar, \
                                  Scalar, Vector, Matrix

from ...utils import InequalityConstraint

import textwrap

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
    Operation(Times(alpha, beta)),
    [], # variants
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(beta, StorageFormat.full)
    ],
    OutputOperand(out, StorageFormat.full), # return value
    cf, # cost function
    "$out = $alpha * $beta",
    [], # Argument objects
    )

# scalar division

alpha = Scalar("alpha")
beta = Scalar("beta")
out = Scalar("out")
cf = lambda d: 1

scalar_division = KernelDescription(
    OperationKV(
        "side",
        {"R": Times(alpha, Inverse(beta)),
         "L": Times(Inverse(beta), alpha)}
        ),
    [], # variants
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(beta, StorageFormat.full)
    ],
    OutputOperand(out, StorageFormat.full), # return value
    cf, # cost function
    "$out = $alpha / $beta",
    [], # Argument objects
    options={KernelOption.no_simplifications}
    )

# scalar inversion

alpha = Scalar("alpha")
out = Scalar("out")
cf = lambda d: 1

scalar_inversion = KernelDescription(
    Operation(Inverse(alpha)),
    [], # variants
    [InputOperand(alpha, StorageFormat.full),
    ],
    OutputOperand(out, StorageFormat.full), # return value
    cf, # cost function
    "$out = 1.0 / $alpha",
    [], # Argument objects
    )

# scalar sum

alpha = Scalar("alpha")
beta = Scalar("beta")
out = Scalar("out")
cf = lambda d: 1

scalar_sum = KernelDescription(
    Operation(Plus(alpha, beta)),
    [], # variants
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(beta, StorageFormat.full)
    ],
    OutputOperand(out, StorageFormat.full), # return value
    cf, # cost function
    "$out = $alpha + $beta",
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
    Operation(Times(Transpose(x),y)),
    [], # variants
    [InputOperand(x, StorageFormat.full),
     InputOperand(y, StorageFormat.full)
    ],
    OutputOperand(out, StorageFormat.full), # return value
    cf, # cost function
    "$out = BLAS.dot($N, $x, 1, $y, 1)",
    [SizeArgument("N", x, "rows")], # Argument objects
    )

# SCAL

x = Vector("x", (n, 1))
alpha = Scalar("alpha")
cf = lambda d: d["N"]

scal = KernelDescription(
    OperationKV(
        None,
        {"L": Times(alpha, x),
         "R": Times(x, alpha)}
    ),
    [], # variants
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(x, StorageFormat.full)
    ],
    OutputOperand(x, StorageFormat.full), # return value
    cf, # cost function
    "scal!($N, $alpha, $x, 1)",
    [SizeArgument("N", x, "rows")], # Argument objects
    options={KernelOption.no_simplifications}
    )

# AXPY

x = Vector("x", (n, 1))
y = Vector("y", (n, 1))
alpha = Scalar("alpha")
cf = lambda d: max(d["N"], d["M"]) # to get correct cost for row and column vectors

axpy = KernelDescription(
    Operation(Plus(Times(alpha, x), y)), # Plus(Times(alpha, Op1(x)), Op2(y))
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
    "axpy!($alpha, $x, $y) # vectors",
    [SizeArgument("N", x, "rows"),
     SizeArgument("M", x, "columns")], # Argument objects
    options={KernelOption.transpose}
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
    Operation(Plus(Times(alpha, x, Transpose(y)), A)),
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
    "ger!($alpha, $x, $y, $A)",
    [SizeArgument("M", x, "rows"),
     SizeArgument("N", y, "rows")], # Argument objects
    )

A = Matrix("A", (m, n))
x = Vector("x", (m, 1))
y = Vector("y", (n, 1))
alpha = Scalar("alpha")
cf = lambda d: 2*d["N"]*d["M"]

ger_alt = KernelDescription(
    Operation(Times(alpha, x, Transpose(y))),
    [
        DefaultValueKV(alpha, [ConstantScalar(1.0)])
    ],
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(x, StorageFormat.full),
     InputOperand(y, StorageFormat.full)
    ],
    OutputOperand(A, StorageFormat.full), # return value
    cf, # cost function
    "$A .= $alpha.*$x.*transpose($y)",
    # textwrap.dedent(
    #     """\
    #     fill!($A, 0.0)
    #     ger!($alpha, $x, $y, $A)\
    #     """
    #     ),
    [SizeArgument("M", x, "rows"),
     SizeArgument("N", y, "rows")], # Argument objects
    )


# SYR (outer product)

A = Matrix("A", (m, n))
x = Vector("x", (m, 1))
alpha = Scalar("alpha")
cf = lambda d: d["M"]**2

syr = KernelDescription(
    Operation(Plus(Times(alpha, x, Transpose(x)), A)),
    [
        DefaultValueKV(alpha, [ConstantScalar(1.0)])
    ],
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(x, StorageFormat.full),
     InputOperand(A, StorageFormat.full)
    ],
    OutputOperand(A, StorageFormat.symmetric_lower_triangular), # return value
    cf, # cost function
    "syr!('L', $alpha, $x, $A)",
    [SizeArgument("M", x, "rows")], # Argument objects
    )

A = Matrix("A", (m, n))
x = Vector("x", (m, 1))
alpha = Scalar("alpha")
cf = lambda d: d["M"]**2

syr_alt = KernelDescription(
    Operation(Times(alpha, x, Transpose(x))),
    [
        DefaultValueKV(alpha, [ConstantScalar(1.0)])
    ],
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(x, StorageFormat.full),
    ],
    OutputOperand(A, StorageFormat.symmetric_lower_triangular), # return value
    cf, # cost function
    # If A is not allocated, use "$A = zeros($type, ($M, $N))\n"
    textwrap.dedent(
        """\
        fill!($A, 0.0)
        syr!('L', $alpha, $x, $A)\
        """
        ),
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
    Operation(Plus(Times(alpha, Op1(A), x), Times(beta, y))),
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
    "gemv!($transA, $alpha, $A, $x, $beta, $y)",
    [SizeArgument("M", A, "rows"),
     SizeArgument("N", A, "columns")], # Argument objects
    options={KernelOption.transpose}
    )


# SYMV

A = Matrix("A", (m, n))
A.set_property(Property.SYMMETRIC)
x = Vector("x", (n, 1))
y = Vector("y", (m, 1))
alpha = Scalar("alpha")
beta = Scalar("beta")
cf = lambda d: 2*d["M"]**2

symv = KernelDescription(
    Operation(Plus(Times(alpha, A, x), Times(beta, y))),
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
    "symv!($uplo, $alpha, $A, $x, $beta, $y)",
    [SizeArgument("M", A, "rows"),
     StorageFormatArgument("uplo", A, {StorageFormat.symmetric_lower_triangular: "L", StorageFormat.symmetric_upper_triangular: "U"}),], # Argument objects
    options={KernelOption.transpose}
    )


# TRMV

"""
TODO: Julia documentation wrong
"""

A = Matrix("A", (n, n))
A.set_property(Property.SQUARE)
x = Vector("x", (n, 1))
cf = lambda d: d["N"]**2

trmv = KernelDescription(
    Operation(Times(Op1(A), x)),
    [
        OperatorKV("transA", {"N": Identity, "T": Transpose}, Op1),
        PropertyKV("uplo", {"U": Property.UPPER_TRIANGULAR, "L": Property.LOWER_TRIANGULAR}, A)
    ],
    [InputOperand(A, StorageFormat.triangular_udiag_opt),
     InputOperand(x, StorageFormat.full),
    ],
    OutputOperand(x, StorageFormat.full), # return value
    cf, # cost function
    "trmv!($uplo, $transA, $diag, $A, $x)",
    [SizeArgument("N", A, "rows"),
     PropertyArgument("diag", A, Property.UNIT_DIAGONAL, ["U", "N"])], # Argument objects
    options={KernelOption.transpose}
    )


# TRSV

A = Matrix("A", (n, n))
A.set_property(Property.SQUARE)
x = Vector("x", (n, 1))
cf = lambda d: d["N"]**2

"""
TODO also use the following?
- trsv(ul, tA, dA, A, b)
Compare performance.
"""

trsv = KernelDescription(
    Operation(Times(Op1(Inverse(A)), x)),
    [
        OperatorKV("transA", {"N": Identity, "T": Transpose}, Op1),
        PropertyKV("uplo", {"U": Property.UPPER_TRIANGULAR, "L": Property.LOWER_TRIANGULAR}, A)
    ],
    [InputOperand(A, StorageFormat.triangular_udiag_opt),
     InputOperand(x, StorageFormat.full),
    ],
    OutputOperand(x, StorageFormat.full), # return value
    cf, # cost function
    "trsv!($uplo, $transA, $diag, $A, $x)",
    [SizeArgument("N", A, "rows"),
     PropertyArgument("diag", A, Property.UNIT_DIAGONAL, ["U", "N"])], # Argument objects
    options={KernelOption.transpose}
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
    Operation(Plus(Times(alpha, Op1(A), Op2(B)), Times(beta, C))),
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
    "gemm!($transA, $transB, $alpha, $A, $B, $beta, $C)",
    [SizeArgument("M", Op1(A), "rows"),
     SizeArgument("N", Op2(B), "columns"),
     SizeArgument("K", Op1(A), "columns")], # Argument objects
    constraints=[InequalityConstraint("A", "C"), InequalityConstraint("B", "C")]
    )


# SYMM

A = Matrix("A", (m, k))
A.set_property(Property.SYMMETRIC)
B = Matrix("B", (k, n))
C = Matrix("C", (m, n))
alpha = Scalar("alpha")
beta = Scalar("beta")
cf = lambda d: 2*d["M"]*d["N"]*d["K"]
"""
The names for the sizes in this cost function do not make a lot of sense.
However, it is not possible to do something more consistent because of the
side argument. A more intuitive way to describe the cost would be n^2*m, where
n is the size of matrix A, and m is the "other" size.
"""

symm = KernelDescription(
    OperationKV(
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
    "symm!($side, $uplo, $alpha, $A, $B, $beta, $C)",
    [SizeArgument("M", B, "rows"),
     SizeArgument("N", B, "columns"),
     SizeArgument("K", A, "rows"),
     StorageFormatArgument("uplo", A, {StorageFormat.symmetric_lower_triangular: "L", StorageFormat.symmetric_upper_triangular: "U"})], # Argument objects
    options={KernelOption.transpose}
    )



# SYRK

A = Matrix("A", (n, k))
C = Matrix("C", (n, n))
C.set_property(Property.SYMMETRIC)
alpha = Scalar("alpha")
beta = Scalar("beta")
cf = lambda d: d["K"]*(d["N"]**2)


syrk = KernelDescription(
    Operation(Plus(Times(alpha, Op1(A), Op1(Transpose(A))), Times(beta, C))),
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
    "syrk!($uplo, $transA, $alpha, $A, $beta, $C)",
    [SizeArgument("N", Op1(A), "rows"),
     SizeArgument("K", Op1(A), "columns"),
     StorageFormatArgument("uplo", C, {StorageFormat.symmetric_lower_triangular: "L", StorageFormat.symmetric_upper_triangular: "U"})], # Argument objects
    )


# SYR2K

A = Matrix("A", (n, k))
B = Matrix("B", (n, k))
C = Matrix("C", (n, n))
C.set_property(Property.SYMMETRIC)
alpha = Scalar("alpha")
beta = Scalar("beta")
cf = lambda d: 2*d["K"]*(d["N"]**2)

syr2k = KernelDescription(
    Operation(Plus(Times(alpha, Op1(A), Op1(Transpose(B))), Times(alpha, Op1(B), Op1(Transpose(A))), Times(beta, C))),
    [
        OperatorKV("trans", {"N": Identity, "T": Transpose}, Op1),
        DefaultValueKV(alpha, [ConstantScalar(1.0)]),
        DefaultValueKV(beta, [ConstantScalar(0.0), ConstantScalar(1.0)])
    ],
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(A, StorageFormat.full),
     InputOperand(B, StorageFormat.full),
     InputOperand(beta, StorageFormat.full),
     InputOperand(C, StorageFormat.symmetric_triangular),
    ],
    OutputOperand(C, StorageFormat.symmetric_triangular_out), # return value
    cf, # cost function
    "syr2k!($uplo, $trans, $alpha, $A, $B, $beta, $C)",
    [SizeArgument("N", Op1(A), "rows"),
     SizeArgument("K", Op1(A), "columns"),
     StorageFormatArgument("uplo", C, {StorageFormat.symmetric_lower_triangular: "L", StorageFormat.symmetric_upper_triangular: "U"})], # Argument objects
    )




# TRMM

A = Matrix("A", (m, m))
A.set_property(Property.SQUARE)
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
    OperationKV(
        "side",
        {"L": Times(alpha, Op1(A), B),
         "R": Times(alpha, B, Op1(A))}
    ),
    [
        OperatorKV("transA", {"N": Identity, "T": Transpose}, Op1),
        DefaultValueKV(alpha, [ConstantScalar(1.0)]),
        PropertyKV("uplo", {"U": Property.UPPER_TRIANGULAR, "L": Property.LOWER_TRIANGULAR}, A)
    ],
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(A, StorageFormat.triangular_udiag_opt),
     InputOperand(B, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    "trmm!($side, $uplo, $transA, $diag, $alpha, $A, $B)",
    [SizeArgument("M", B, "rows"),
     SizeArgument("N", B, "columns"),
     SizeArgument("K", A, "rows"),
     PropertyArgument("diag", A, Property.UNIT_DIAGONAL, ["U", "N"])], # Argument objects
    options={KernelOption.transpose}
    )


# TRSM

A = Matrix("A", (m, m))
A.set_property(Property.SQUARE)
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
    OperationKV(
        "side",
        {"L": Times(alpha, Op1(Inverse(A)), B),
         "R": Times(alpha, B, Op1(Inverse(A)))}
    ),
    [
        OperatorKV("transA", {"N": Identity, "T": Transpose}, Op1),
        DefaultValueKV(alpha, [ConstantScalar(1.0)]),
        PropertyKV("uplo", {"U": Property.UPPER_TRIANGULAR, "L": Property.LOWER_TRIANGULAR}, A)
    ],
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(A, StorageFormat.triangular_udiag_opt),
     InputOperand(B, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    "trsm!($side, $uplo, $transA, $diag, $alpha, $A, $B)",
    [SizeArgument("M", B, "rows"),
     SizeArgument("N", B, "columns"),
     SizeArgument("K", A, "rows"),
     PropertyArgument("diag", A, Property.UNIT_DIAGONAL, ["U", "N"])], # Argument objects
    options={KernelOption.transpose}
    )

###############
# LAPACK

# SCAL (scalar * matrix) (equivalent to LASCL in LAPACK)

X = Matrix("X", (m, n))
alpha = Scalar("alpha")
cf = lambda d: d["M"]*d["N"]

lascl = KernelDescription(
    OperationKV(
        None,
        {"L": Times(alpha, X),
         "R": Times(X, alpha)}
    ),
    [], # variants
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(X, StorageFormat.full),
    ],
    OutputOperand(X, StorageFormat.full), # return value
    cf, # cost function
    "scal!($MN, $alpha, $X, 1)",
    [SizeArgument("M", X, "rows"),
     SizeArgument("N", X, "columns"),
     SizeArgument("MN", X, "entries")], # Argument objects
    options={KernelOption.no_simplifications}
    )


# getri (matrix inversion)

A = Matrix("A", (n, n))
B = Matrix("B", (n, n))
A.set_property(Property.SQUARE)
cf = lambda d: 2*d["N"]**3

getri = KernelDescription(
    Operation(Inverse(A)),
    [], # variants
    [InputOperand(A, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    "$B = inv($A)",
    [SizeArgument("N", A, "columns")], # Argument objects
    options={KernelOption.transpose}
    )

# trtri (triangular matrix inversion)

A = Matrix("A", (n, n))
A.set_property(Property.SQUARE)
A.set_property(Property.TRIANGULAR)
cf = lambda d: (d["N"]**3)/3

trtri = KernelDescription(
    Operation(Inverse(A)),
    [
        PropertyKV("uplo", {"U": Property.UPPER_TRIANGULAR, "L": Property.LOWER_TRIANGULAR}, A)
    ], # variants
    [InputOperand(A, StorageFormat.triangular_udiag_opt),
    ],
    OutputOperand(A, StorageFormat.as_overwritten), # return value
    cf, # cost function
    "LAPACK.trtri!($uplo, $diag, $A)",
    [SizeArgument("N", A, "rows"),
     PropertyArgument("diag", A, Property.UNIT_DIAGONAL, ["U", "N"])], # Argument objects
    options={KernelOption.transpose}
    )


# symmetric matrix inversion)

A = Matrix("A", (n, n))
A.set_property(Property.SQUARE)
A.set_property(Property.SYMMETRIC)
cf = lambda d: (d["N"]**3)/3

syinv = KernelDescription(
    Operation(Inverse(A)),
    [], # variants
    [InputOperand(A, StorageFormat.symmetric_triangular),
    ],
    OutputOperand(A, StorageFormat.symmetric_triangular_out), # return value
    cf, # cost function
    textwrap.dedent(
        """\
        ($A, ipiv, info) = LAPACK.sytrf!($uplo, $A)
        LAPACK.sytri!($uplo, $A, ipiv)\
        """
        ),
    [SizeArgument("N", A, "rows"),
     StorageFormatArgument("uplo", A, {StorageFormat.symmetric_lower_triangular: "L", StorageFormat.symmetric_upper_triangular: "U"})], # Argument objects
    options={KernelOption.transpose}
    )


# POSV

A = Matrix("A", (m, m))
A.set_property(Property.SQUARE)
A.set_property(Property.SPD)
B = Matrix("B", (m, n))
cf = lambda d: (d["M"]**3)/3 + 2*(d["M"]**2)*d["N"]

"""
TODO problem: both A and B are overwritten, but it's not possible to express that here
alternative
LAPACK.potrf!('L', A)
LAPACK.potrs!('L', A, B)
"""

posv = KernelDescription(
    Operation(Times(Inverse(A), B)),
    [],  # variants
    [InputOperand(A, StorageFormat.symmetric_triangular),
     InputOperand(B, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    "LAPACK.posv!('L', $A, $B)",
    [SizeArgument("M", B, "rows"),
     SizeArgument("N", B, "columns")], # Argument objects
    options={KernelOption.transpose}
    )


# POSV right

A = Matrix("A", (m, m))
A.set_property(Property.SQUARE)
A.set_property(Property.SPD)
B = Matrix("B", (n, m))
cf = lambda d: (d["M"]**3)/3 + 2*(d["M"]**2)*d["N"]

"""
TODO problem: both A and B are overwritten, but it's not possible to express that here
alternative
LAPACK.potrf!('L', A)
LAPACK.potrs!('L', A, B)
"""

posvr = KernelDescription(
    Operation(Times(B, Inverse(A))),
    [],  # variants
    [InputOperand(A, StorageFormat.symmetric_triangular),
     InputOperand(B, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    textwrap.dedent(
        """\
        $B = $B'
        LAPACK.posv!('L', $A, $B)
        $B = $B'\
        """
        ),
    [SizeArgument("M", B, "columns"),
     SizeArgument("N", B, "rows")], # Argument objects
    options={KernelOption.transpose}
    )


# sysv

A = Matrix("A", (m, m))
A.set_property(Property.SQUARE)
A.set_property(Property.SYMMETRIC)
B = Matrix("B", (m, n))
cf = lambda d: (d["M"]**3)/3 + 2*(d["M"]**2)*d["N"]

"""
TODO problem: both A and B are overwritten, but it's not possible to express that here
alternative. As a preliminary solution, we copy A.

(A, ipiv) = LAPACK.sytrf!('L', A)
LAPACK.sytrs!('L', A, ipiv, B)
TODO For whatever reason, sytrs is very slow. Investigate.
"""

sysv = KernelDescription(
    Operation(Times(Inverse(A), B)),
    [],  # variants
    [InputOperand(A, StorageFormat.symmetric_triangular),
     InputOperand(B, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    textwrap.dedent(
        """\
        tmp = Array{Float64}(undef, $M, $M)
        blascopy!($M*$M, $A, 1, tmp, 1)
        LAPACK.sysv!('L', tmp, $B)\
        """
        ),
    [SizeArgument("M", B, "rows"),
     SizeArgument("N", B, "columns")], # Argument objects
    options={KernelOption.transpose}
    )


# sysv right

A = Matrix("A", (m, m))
A.set_property(Property.SQUARE)
A.set_property(Property.SYMMETRIC)
B = Matrix("B", (n, m))
cf = lambda d: (d["M"]**3)/3 + 2*(d["M"]**2)*d["N"]

"""
TODO problem: both A and B are overwritten, but it's not possible to express that here
alternative
(A, ipiv) = LAPACK.sytrf!('L', A)
LAPACK.sytrs!('L', A, ipiv, B)
TODO For whatever reason, sytrs is very slow. Investigate.
"""

sysvr = KernelDescription(
    Operation(Times(B, Inverse(A))),
    [],  # variants
    [InputOperand(A, StorageFormat.symmetric_triangular),
     InputOperand(B, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    textwrap.dedent(
        """\
        $B = transpose($B)
        tmp = Array{Float64}(undef, $M, $M)
        blascopy!($M*$M, $A, 1, tmp, 1)
        LAPACK.sysv!('L', tmp, $B)
        $B = transpose($B)\
        """
        ),
    [SizeArgument("M", B, "columns"),
     SizeArgument("N", B, "rows")], # Argument objects
    options={KernelOption.transpose}
    )

# gesv

A = Matrix("A", (m, m))
A.set_property(Property.SQUARE)
B = Matrix("B", (m, n))
cf = lambda d: 2*(d["M"]**3)/3 + 2*(d["M"]**2)*d["N"]

"""
TODO problem: both A and B are overwritten, but it's not possible to express that here
alternative
LAPACK.gesv!($A, $B)
"""

gesv = KernelDescription(
    Operation(Times(Op1(Inverse(A)), B)),
    [
        OperatorKV("transA", {"N": Identity, "T": Transpose}, Op1),
    ],  # variants
    [InputOperand(A, StorageFormat.full),
     InputOperand(B, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    "($A, ipiv, info) = LAPACK.getrf!($A)\nLAPACK.getrs!($transA, $A, ipiv, $B)",
    [SizeArgument("M", B, "rows"),
     SizeArgument("N", B, "columns")], # Argument objects
    options={KernelOption.transpose}
    )


# gesv right

A = Matrix("A", (m, m))
A.set_property(Property.SQUARE)
B = Matrix("B", (n, m))
cf = lambda d: 2*(d["M"]**3)/3 + 2*(d["M"]**2)*d["N"]

"""
TODO problem: both A and B are overwritten, but it's not possible to express that here
alternative
LAPACK.gesv!($A, $B)
"""

gesvr = KernelDescription(
    Operation(Times(B, Inverse(A))),
    [],  # variants
    [InputOperand(A, StorageFormat.full),
     InputOperand(B, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    "$B = $B/lufact!($A)",
    [SizeArgument("M", B, "columns"),
     SizeArgument("N", B, "rows")], # Argument objects
    options={KernelOption.transpose}
    )

# gesv right transpose

A = Matrix("A", (m, m))
A.set_property(Property.SQUARE)
B = Matrix("B", (n, m))
cf = lambda d: 2*(d["M"]**3)/3 + 2*(d["M"]**2)*d["N"]

"""
TODO problem: both A and B are overwritten, but it's not possible to express that here
alternative
LAPACK.gesv!($A, $B)
"""

gesvrt = KernelDescription(
    Operation(Times(B, InverseTranspose(A))),
    [],  # variants
    [InputOperand(A, StorageFormat.full),
     InputOperand(B, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    "$B = $B/lufact!($A')",
    [SizeArgument("M", B, "columns"),
     SizeArgument("N", B, "rows")], # Argument objects
    options={KernelOption.transpose}
    )


#################
#################

# POSV vector

A = Matrix("A", (m, m))
A.set_property(Property.SQUARE)
A.set_property(Property.SPD)
x = Vector("x", (m, 1))
cf = lambda d: (d["M"]**3)/3 + 2*(d["M"]**2)

"""
TODO problem: both A and B are overwritten, but it's not possible to express that here
alternative
LAPACK.potrf!('L', A)
LAPACK.potrs!('L', A, B)
"""

posv_vec = KernelDescription(
    Operation(Times(Inverse(A), x)),
    [],  # variants
    [InputOperand(A, StorageFormat.symmetric_triangular),
     InputOperand(x, StorageFormat.full),
    ],
    OutputOperand(x, StorageFormat.full), # return value
    cf, # cost function
    "LAPACK.posv!('L', $A, $x)",
    [SizeArgument("M", A, "rows")], # Argument objects
    options={KernelOption.transpose}
    )


# # POSV right vector
# # probably not necessary

# A = Matrix("A", (m, m))
# A.set_property(Property.SQUARE)
# A.set_property(Property.SPD)
# x = Vector("x", (m, 1))
# cf = lambda d: (d["M"]**3)/3 + 2*(d["M"]**2)*d["N"]

# """
# TODO problem: both A and B are overwritten, but it's not possible to express that here
# alternative
# LAPACK.potrf!('L', A)
# LAPACK.potrs!('L', A, B)
# """

# posvr = KernelDescription(
#     OperationKV(
#         None,
#         {None: Times(Transpose(x), Inverse(A))}
#     ),
#     [],  # variants
#     [InputOperand(A, StorageFormat.symmetric_triangular),
#      InputOperand(B, StorageFormat.full),
#     ],
#     OutputOperand(B, StorageFormat.full), # return value
#     cf, # cost function
#
#     "$B = $B'\nLAPACK.posv!('L', $A, $B)\n$B = $B'",
#
#     [SizeArgument("M", B, "columns"),
#      SizeArgument("N", B, "rows")], # Argument objects
#     options={KernelOption.transpose}
#     )


# sysv vector

A = Matrix("A", (m, m))
A.set_property(Property.SQUARE)
A.set_property(Property.SYMMETRIC)
x = Vector("x", (m, 1))
cf = lambda d: (d["M"]**3)/3 + 2*(d["M"]**2)

"""
TODO problem: both A and B are overwritten, but it's not possible to express that here
alternative
(A, ipiv) = LAPACK.sytrf!('L', A)
LAPACK.sytrs!('L', A, ipiv, B)
TODO For whatever reason, sytrs is very slow. Investigate.
"""

sysv_vec = KernelDescription(
    Operation(Times(Inverse(A), x)),
    [],  # variants
    [InputOperand(A, StorageFormat.symmetric_triangular),
     InputOperand(x, StorageFormat.full),
    ],
    OutputOperand(x, StorageFormat.full), # return value
    cf, # cost function
    textwrap.dedent(
        """\
        tmp = Array{Float64}(undef, $M, $M)
        blascopy!($M*$M, $A, 1, tmp, 1)
        LAPACK.sysv!('L', tmp, $x)\
        """
        ),
    [SizeArgument("M", A, "rows")], # Argument objects
    options={KernelOption.transpose}
    )


# # sysv right vector
# # probably not necessary

# A = Matrix("A", (m, m))
# A.set_property(Property.SQUARE)
# A.set_property(Property.SYMMETRIC)
# x = Vector("x", (m, 1))
# cf = lambda d: (d["M"]**3)/3 + 2*(d["M"]**2)*d["N"]

# """
# TODO problem: both A and B are overwritten, but it's not possible to express that here
# alternative
# (A, ipiv) = LAPACK.sytrf!('L', A)
# LAPACK.sytrs!('L', A, ipiv, B)
# TODO For whatever reason, sytrs is very slow. Investigate.
# """

# sysvr = KernelDescription(
#     OperationKV(
#         None,
#         {None: Times(Transpose(x), Inverse(A))}
#     ),
#     [],  # variants
#     [InputOperand(A, StorageFormat.symmetric_triangular),
#      InputOperand(B, StorageFormat.full),
#     ],
#     OutputOperand(B, StorageFormat.full), # return value
#     cf, # cost function
#
#     "$B = $B'\nLAPACK.sysv!('L', $A, $B)\n$B = $B'",
#
#     [SizeArgument("M", B, "columns"),
#      SizeArgument("N", B, "rows")], # Argument objects
#     options={KernelOption.transpose}
#     )

# gesv vector

A = Matrix("A", (m, m))
A.set_property(Property.SQUARE)
x = Vector("x", (m, 1))
cf = lambda d: 2*(d["M"]**3)/3 + 2*(d["M"]**2)

"""
TODO problem: both A and B are overwritten, but it's not possible to express that here
alternative
LAPACK.gesv!($A, $B)
"""

gesv_vec = KernelDescription(
    Operation(Times(Op1(Inverse(A)), x)),
    [
        OperatorKV("transA", {"N": Identity, "T": Transpose}, Op1),
    ],  # variants
    [InputOperand(A, StorageFormat.full),
     InputOperand(x, StorageFormat.full),
    ],
    OutputOperand(x, StorageFormat.full), # return value
    cf, # cost function
    textwrap.dedent(
        """\
        ($A, ipiv, info) = LAPACK.getrf!($A)
        LAPACK.getrs!($transA, $A, ipiv, $x)\
        """
        ),
    [SizeArgument("M", A, "rows")], # Argument objects
    options={KernelOption.transpose}
    )


# # gesv right vector
# # probably not necessary

# A = Matrix("A", (m, m))
# A.set_property(Property.SQUARE)
# x = Vector("x", (m, 1))
# cf = lambda d: 2*(d["M"]**3)/3 + 2*(d["M"]**2)*d["N"]

# """
# TODO problem: both A and B are overwritten, but it's not possible to express that here
# alternative
# LAPACK.gesv!($A, $B)
# """

# gesvr = KernelDescription(
#     OperationKV(
#         None,
#         {None: Times(B, Inverse(A))}
#     ),
#     [],  # variants
#     [InputOperand(A, StorageFormat.full),
#      InputOperand(B, StorageFormat.full),
#     ],
#     OutputOperand(B, StorageFormat.full), # return value
#     cf, # cost function
#
#     "$B = $B/lufact!($A)",
#
#     [SizeArgument("M", B, "columns"),
#      SizeArgument("N", B, "rows")], # Argument objects
#     options={KernelOption.transpose}
#     )

# # gesv right transpose vector
# # probably not necessary

# A = Matrix("A", (m, m))
# A.set_property(Property.SQUARE)
# x = Vector("x", (m, 1))
# cf = lambda d: 2*(d["M"]**3)/3 + 2*(d["M"]**2)*d["N"]

# """
# TODO problem: both A and B are overwritten, but it's not possible to express that here
# alternative
# LAPACK.gesv!($A, $B)
# """

# gesvrt = KernelDescription(
#     OperationKV(
#         None,
#         {None: Times(B, InverseTranspose(A))}
#     ),
#     [],  # variants
#     [InputOperand(A, StorageFormat.full),
#      InputOperand(B, StorageFormat.full),
#     ],
#     OutputOperand(B, StorageFormat.full), # return value
#     cf, # cost function
#
#     "$B = $B/lufact!($A')",
#
#     [SizeArgument("M", B, "columns"),
#      SizeArgument("N", B, "rows")], # Argument objects
#     options={KernelOption.transpose}
#     )

#################
#################

# inv(scalar) * vector

x = Vector("x", (n, 1))
alpha = Scalar("alpha")
cf = lambda d: d["N"]

invscal = KernelDescription(
    OperationKV(
        None,
        {"L": Times(Inverse(alpha), x),
         "R": Times(x, Inverse(alpha))}
    ),
    [], # variants
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(x, StorageFormat.full)
    ],
    OutputOperand(x, StorageFormat.full), # return value
    cf, # cost function
    "$x ./= $alpha",
    [SizeArgument("N", x, "rows")], # Argument objects
    options={KernelOption.no_simplifications}
    )

# inv(scalar) * matrix

X = Matrix("X", (m, n))
alpha = Scalar("alpha")
cf = lambda d: d["M"]*d["N"]

invlascl = KernelDescription(
    OperationKV(
        None,
        {"L": Times(Inverse(alpha), X),
         "R": Times(X, Inverse(alpha))}
    ),
    [], # variants
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(X, StorageFormat.full),
    ],
    OutputOperand(X, StorageFormat.full), # return value
    cf, # cost function
    "$X ./= $alpha",
    [SizeArgument("M", X, "rows"),
     SizeArgument("N", X, "columns")], # Argument objects
    options={KernelOption.no_simplifications}
    )

# diaginv (diagonal matrix inversion)

A = Matrix("A", (n, n))
A.set_property(Property.SQUARE)
A.set_property(Property.DIAGONAL)
cf = lambda d: d["N"]

diaginv = KernelDescription(
    Operation(Inverse(A)),
    [], # variants
    [InputOperand(A, StorageFormat.diagonal_vector),
    ],
    OutputOperand(A, StorageFormat.diagonal_vector), # return value
    cf, # cost function
    "$A = 1 ./$A",
    [SizeArgument("N", A, "rows")], # Argument objects
    options={KernelOption.transpose}
    )

# matrix transposition

A = Matrix("A", (m, n))
B = Matrix("B", (n, m))
cf = lambda d: 1

transpose = KernelDescription(
    Operation(Transpose(A)),
    [], # variants
    [InputOperand(A, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    # $B = Array{$type}($N, $M)
    """transpose!($B, $A)""",
    [SizeArgument("M", A, "rows"),
     SizeArgument("N", A, "columns")], # Argument objects
    )

# vector transposition (totally fake) (this is kind of tricky. Since it doesn't actually do anything, overwriting is never a problem. Solution: Just use different symbol for output?)

x = Vector("x", (n, 1))
cf = lambda d: 0

transpose_vector = KernelDescription(
    Operation(Transpose(x)),
    [], # variants
    [InputOperand(x, StorageFormat.full),
    ],
    OutputOperand(x, StorageFormat.full), # return value
    cf, # cost function
    "# Transposing vector $x (no operation);",
    [], # Argument objects
    )

# matrix sum

alpha = Scalar("alpha")
A = Matrix("A", (m, n))
B = Matrix("B", (m, n))
cf = lambda d: 2*d["N"]*d["M"]

matrix_sum = KernelDescription(
    Operation(Plus(Times(alpha, A), B)),
    [
        DefaultValueKV(alpha, [ConstantScalar(1.0)])
    ],
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(A, StorageFormat.full),
     InputOperand(B, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    "axpy!($alpha, $A, $B) # matrices",
    [SizeArgument("N", A, "columns"),
     SizeArgument("M", A, "rows")], # Argument objects
    options={KernelOption.transpose}
    )

# A + B^T

A = Matrix("A", (m, n))
B = Matrix("B", (n, m))
cf = lambda d: d["N"]*d["M"]

matrix_sum_transpose = KernelDescription(
    Operation(Plus(A, Transpose(B))),
    [],
    [InputOperand(A, StorageFormat.full),
     InputOperand(B, StorageFormat.full),
    ],
    OutputOperand(A, StorageFormat.full), # return value
    cf, # cost function
    "$A .+= transpose($B)",
    [SizeArgument("N", A, "columns"),
     SizeArgument("M", A, "rows")], # Argument objects
    )

# scalar * diagonal

X = Matrix("X", (m, n))
X.set_property(Property.DIAGONAL)
alpha = Scalar("alpha")
cf = lambda d: min(d["M"], d["N"])

diagscal = KernelDescription(
    OperationKV(
        None,
        {"L": Times(alpha, X),
         "R": Times(X, alpha)}
    ),
    [], # variants
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(X, StorageFormat.diagonal_vector),
    ],
    OutputOperand(X, StorageFormat.diagonal_vector), # return value
    cf, # cost function
    "scal!(min($M, $N), $alpha, $X, 1)",
    [SizeArgument("M", X, "rows"),
     SizeArgument("N", X, "columns")], # Argument objects
    options={KernelOption.no_simplifications}
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
A.set_property(Property.DIAGONAL)
B = Matrix("B", (k, n))
B.set_property(Property.DIAGONAL)
C = Matrix("C", (m, n))
cf = lambda d: min(d["M"], d["K"], d["N"])

diagdiagmul = KernelDescription(
    Operation(Times(Op1(A), Op2(B))),
    [
        OperatorKV("", {"N": Identity, "T": Transpose}, Op1),
        OperatorKV("", {"N": Identity, "T": Transpose}, Op2),
    ], # variants
    [InputOperand(A, StorageFormat.diagonal_vector),
     InputOperand(B, StorageFormat.diagonal_vector),
    ],
    OutputOperand(C, StorageFormat.diagonal_vector), # return value
    cf, # cost function
    textwrap.dedent(
        """\
        for i = 1:min(length($A), length($B));
            $C[i] = $A[i] * $B[i];
        end;
        if max(length($A), length($B)) < length($C)
            $C[min(length($A), length($B))+1:end] = 0;
        end;\
        """
        ),
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
A.set_property(Property.DIAGONAL)
A.set_property(Property.SQUARE)
B = Matrix("B", (m, n))
B.set_property(Property.DIAGONAL)
cf = lambda d: min(d["M"], d["N"])

diagdiagsolve = KernelDescription(
    OperationKV(
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
    "$B ./= $A;",
    [SizeArgument("M", B, "rows"),
     SizeArgument("N", B, "columns")], # Argument objects
    options={KernelOption.transpose}
    )


# full * inv(diag)

# scaling columns

"""
TODO: what about transB?
"""

A = Matrix("A", (m, m))
A.set_property(Property.DIAGONAL)
A.set_property(Property.SQUARE)
B = Matrix("B", (m, n))
cf = lambda d: d["M"]*d["N"]

diagsmr = KernelDescription(
    Operation(Times(B, Inverse(A))),
    [], # variants
    [InputOperand(A, StorageFormat.diagonal_vector),
     InputOperand(B, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    textwrap.dedent(
        """\
        @views for i = 1:size($B, 2);
            $B[:,i] ./= $A[i];
        end;\
        """
        ),
    [SizeArgument("M", B, "rows"),
     SizeArgument("N", B, "columns")], # Argument objects
    options={KernelOption.transpose}
    )


# inv(diag) * full

# scaling rows

"""
TODO: what about transB?
"""

A = Matrix("A", (m, m))
A.set_property(Property.DIAGONAL)
A.set_property(Property.SQUARE)
B = Matrix("B", (m, n))
cf = lambda d: d["M"]*d["N"]

diagsml = KernelDescription(
    Operation(Times(Inverse(A), B)),
    [], # variants
    [InputOperand(A, StorageFormat.diagonal_vector),
     InputOperand(B, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    "$B ./= $A;",
    [SizeArgument("M", B, "rows"),
     SizeArgument("N", B, "columns")], # Argument objects
    options={KernelOption.transpose}
    )


# inv(diag) * vector

"""
TODO: does the transpose variant get simplified correctly? (i.e. A is symmetric)
"""

A = Matrix("A", (m, m))
A.set_property(Property.DIAGONAL)
A.set_property(Property.SQUARE)
x = Vector("x", (m, 1))
cf = lambda d: d["M"]

diagsv = KernelDescription(
    Operation(Times(Inverse(A), x)),
    [], # variants
    [InputOperand(A, StorageFormat.diagonal_vector),
     InputOperand(x, StorageFormat.full),
    ],
    OutputOperand(x, StorageFormat.full), # return value
    cf, # cost function
    "$x ./= $A",
    [SizeArgument("M", A, "rows")], # Argument objects
    options={KernelOption.transpose}
    )


# full * diag (square)

# scaling columns

"""
TODO: what about transB?

scale!(A, b) is deprecated
"""

A = Matrix("A", (m, m))
A.set_property(Property.DIAGONAL)
"""Symmetric used to be square. While square might be more intuitive, it results
in a property set that is not in canonical form. With symmetric, the set is in
canonical form.
"""
A.set_property(Property.SYMMETRIC)
B = Matrix("B", (m, n))
cf = lambda d: d["M"]*d["N"]

diagmmr = KernelDescription(
    Operation(Times(B, A)),
    [], # variants
    [InputOperand(A, StorageFormat.diagonal_vector),
     InputOperand(B, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    textwrap.dedent(
        """\
        for i = 1:size($B, 2);
            view($B, :, i)[:] .*= $A[i];
        end;\
        """
        ),
    [SizeArgument("M", B, "rows"),
     SizeArgument("N", B, "columns")], # Argument objects
    options={KernelOption.transpose}
    )


# diag (square) * full

# scaling rows

"""
TODO: what about transB?

scale!(b, A) is deprecated
"""

A = Matrix("A", (m, m))
A.set_property(Property.DIAGONAL)
"""Symmetric used to be square. While square might be more intuitive, it results
in a property set that is not in canonical form. With symmetric, the set is in
canonical form.
"""
A.set_property(Property.SYMMETRIC)
B = Matrix("B", (m, n))
cf = lambda d: d["M"]*d["N"]

diagmml = KernelDescription(
    Operation(Times(A, B)),
    [], # variants
    [InputOperand(A, StorageFormat.diagonal_vector),
     InputOperand(B, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    textwrap.dedent(
        """\
        for i = 1:size($B, 2);
            view($B, :, i)[:] .*= $A;
        end;\
        """
        ), # this has better spacial locality than view($B, i, :)[:] .*= $A[i];
    [SizeArgument("M", B, "rows"),
     SizeArgument("N", B, "columns")], # Argument objects
    options={KernelOption.transpose}
    )


# diag (square) * vector

"""
TODO: does the transpose variant get simplified correctly? (i.e. A is symmetric)
"""

A = Matrix("A", (m, m))
A.set_property(Property.DIAGONAL)
"""Symmetric used to be square. While square might be more intuitive, it results
in a property set that is not in canonical form. With symmetric, the set is in
canonical form.
"""
A.set_property(Property.SYMMETRIC)
x = Vector("x", (m, 1))
cf = lambda d: d["M"]

diagmv = KernelDescription(
    Operation(Times(A, x)),
    [], # variants
    [InputOperand(A, StorageFormat.diagonal_vector),
     InputOperand(x, StorageFormat.full),
    ],
    OutputOperand(x, StorageFormat.full), # return value
    cf, # cost function
    "$x .*= $A",
    [SizeArgument("M", A, "rows")], # Argument objects
    options={KernelOption.transpose}
    )

# diag + diag

A = Matrix("A", (m, n))
A.set_property(Property.DIAGONAL)
B = Matrix("B", (m, n))
B.set_property(Property.DIAGONAL)
alpha = Scalar("alpha")
cf = lambda d: min(d["M"], d["N"])

diagdiagadd = KernelDescription(
    Operation(Plus(Times(alpha, A), B)),
    [
        DefaultValueKV(alpha, [ConstantScalar(1.0)]),
    ],
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(A, StorageFormat.diagonal_vector),
     InputOperand(B, StorageFormat.diagonal_vector)
    ],
    OutputOperand(B, StorageFormat.diagonal_vector), # return value
    cf, # cost function
    "axpy!($alpha, $A, $B) # diagonal matrices",
    [SizeArgument("N", A, "rows"),
     SizeArgument("M", A, "columns")], # Argument objects
    )


# diag + full

A = Matrix("A", (m, n))
A.set_property(Property.DIAGONAL)
B = Matrix("B", (m, n))
alpha = Scalar("alpha")
cf = lambda d: min(d["M"], d["N"])

diagfulladd = KernelDescription(
    Operation(Plus(Times(alpha, A), B)),
    [
        DefaultValueKV(alpha, [ConstantScalar(1.0)]),
    ],
    [InputOperand(alpha, StorageFormat.full),
     InputOperand(A, StorageFormat.diagonal_vector),
     InputOperand(B, StorageFormat.explicit_diagonal)
    ],
    OutputOperand(B, StorageFormat.as_overwritten), # return value
    cf, # cost function
    textwrap.dedent(
        """\
        d = $A;
        for i=1:length(d)
            $B[i, i] += $alpha*d[i];
        end;\
        """
        ), # the variable d is used here to make sure that the identity matrix is no created inside the loop
    [SizeArgument("N", A, "rows"),
     SizeArgument("M", A, "columns")], # Argument objects
    options={KernelOption.transpose}
    )


# permutation * matrix

P = Matrix("P", (n, n))
P.set_property(Property.PERMUTATION)
A = Matrix("A", (n, m))
B = Matrix("B", (n, m))
cf = lambda d: d["N"]*d["M"]

pmm = KernelDescription(
    Operation(Times(P, A)),
    [], # variants
    [InputOperand(P, StorageFormat.permutation_vector),
     InputOperand(A, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    "$B = $A[$P,:]",
    [SizeArgument("N", A, "rows"),
     SizeArgument("M", A, "columns")], # Argument objects
    options={KernelOption.transpose}
    )


# transpose(permutation) * matrix

P = Matrix("P", (n, n))
P.set_property(Property.PERMUTATION)
A = Matrix("A", (n, m))
B = Matrix("B", (n, m))
cf = lambda d: d["N"]*d["M"]

ptmm = KernelDescription(
    Operation(Times(Transpose(P), A)),
    [], # variants
    [InputOperand(P, StorageFormat.permutation_vector),
     InputOperand(A, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    "$B = $A[invperm($P),:]",
    [SizeArgument("N", A, "rows"),
     SizeArgument("M", A, "columns")], # Argument objects
    options={KernelOption.transpose}
    )


# matrix * permutation

A = Matrix("A", (n, m))
P = Matrix("P", (m, m))
P.set_property(Property.PERMUTATION)
B = Matrix("B", (n, m))
cf = lambda d: d["N"]*d["M"]

mpm = KernelDescription(
    Operation(Times(A, P)),
    [], # variants
    [InputOperand(P, StorageFormat.permutation_vector),
     InputOperand(A, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    "$B = $A[:,invperm($P)]",
    [SizeArgument("N", A, "rows"),
     SizeArgument("M", A, "columns")], # Argument objects
    options={KernelOption.transpose}
    )


# matrix * transpose(permutation)

A = Matrix("A", (n, m))
P = Matrix("P", (m, m))
P.set_property(Property.PERMUTATION)
B = Matrix("B", (n, m))
cf = lambda d: d["N"]*d["M"]

mptm = KernelDescription(
    Operation(Times(A, Transpose(P))),
    [], # variants
    [InputOperand(P, StorageFormat.permutation_vector),
     InputOperand(A, StorageFormat.full),
    ],
    OutputOperand(B, StorageFormat.full), # return value
    cf, # cost function
    "$B = $A[:,$P]",
    [SizeArgument("N", A, "rows"),
     SizeArgument("M", A, "columns")], # Argument objects
    options={KernelOption.transpose}
    )

# permutation * vector

P = Matrix("P", (n, n))
P.set_property(Property.PERMUTATION)
x = Vector("x", (n, 1))
y = Vector("y", (n, 1))
cf = lambda d: d["N"]

pvm = KernelDescription(
    Operation(Times(P, x)),
    [], # variants
    [InputOperand(P, StorageFormat.permutation_vector),
     InputOperand(x, StorageFormat.full),
    ],
    OutputOperand(y, StorageFormat.full), # return value
    cf, # cost function
    "$y = $x[$P]",
    [SizeArgument("N", P, "rows")], # Argument objects
    options={KernelOption.transpose}
    )


# transpose(permutation) * vector

P = Matrix("P", (n, n))
P.set_property(Property.PERMUTATION)
x = Vector("x", (n, 1))
y = Vector("y", (n, 1))
cf = lambda d: d["N"]

ptvm = KernelDescription(
    Operation(Times(Transpose(P), x)),
    [], # variants
    [InputOperand(P, StorageFormat.permutation_vector),
     InputOperand(x, StorageFormat.full),
    ],
    OutputOperand(y, StorageFormat.full), # return value
    cf, # cost function
    "$y = $x[invperm($P)]",
    [SizeArgument("N", P, "rows")], # Argument objects
    options={KernelOption.transpose}
    )

# permutation * permutation

P = Matrix("P", (n, n))
P.set_property(Property.PERMUTATION)
Q = Matrix("Q", (n, n))
Q.set_property(Property.PERMUTATION)
X = Matrix("X", (n, n))
cf = lambda d: 0

ppm = KernelDescription(
    Operation(Times(P, Q)),
    [], # variants
    [InputOperand(P, StorageFormat.permutation_vector),
     InputOperand(Q, StorageFormat.permutation_vector),
    ],
    OutputOperand(X, StorageFormat.permutation_vector), # return value
    cf, # cost function
    "$X = $Q[$P]",
    [SizeArgument("N", P, "rows")], # Argument objects
    options={KernelOption.transpose}
    )


# transpose(permutation) * permutation

P = Matrix("P", (n, n))
P.set_property(Property.PERMUTATION)
Q = Matrix("Q", (n, n))
Q.set_property(Property.PERMUTATION)
X = Matrix("X", (n, n))
cf = lambda d: 0

ptpm = KernelDescription(
    Operation(Times(Transpose(P), Q)),
    [], # variants
    [InputOperand(P, StorageFormat.permutation_vector),
     InputOperand(Q, StorageFormat.permutation_vector),
    ],
    OutputOperand(X, StorageFormat.permutation_vector), # return value
    cf, # cost function
    "$X = $Q[invperm($P)]",
    [SizeArgument("N", P, "rows")], # Argument objects
    )


# permutation * transpose(permutation)

P = Matrix("P", (n, n))
P.set_property(Property.PERMUTATION)
Q = Matrix("Q", (n, n))
Q.set_property(Property.PERMUTATION)
X = Matrix("X", (n, n))
cf = lambda d: 0

pptm = KernelDescription(
    Operation(Times(P, Transpose(Q))),
    [], # variants
    [InputOperand(P, StorageFormat.permutation_vector),
     InputOperand(Q, StorageFormat.permutation_vector),
    ],
    OutputOperand(X, StorageFormat.permutation_vector), # return value
    cf, # cost function
    "$X = invperm($Q)[$P]",
    [SizeArgument("N", P, "rows")], # Argument objects
    )


# transpose(permutation)

P = Matrix("P", (n, n))
P.set_property(Property.PERMUTATION)
Q = Matrix("Q", (n, n))
cf = lambda d: 0

transpose_perm = KernelDescription(
    Operation(Transpose(P)),
    [], # variants
    [InputOperand(P, StorageFormat.permutation_vector),
    ],
    OutputOperand(Q, StorageFormat.permutation_vector), # return value
    cf, # cost function
    "$Q = invperm($P)",
    [SizeArgument("N", P, "rows")], # Argument objects
    )

# what about diagonal A? could be almost the same as vector

