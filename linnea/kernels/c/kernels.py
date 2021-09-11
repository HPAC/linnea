# from clak.kernels.utils.kernels import KernelDescription, \
#                                           Op1, Op2, \
#                                           OperationKV, PropertyKV, DefaultValueKV, OperatorKV, \
#                                           KernelType
# from clak.kernels.utils.general import SizeArgument, PropertyArgument, StrideArgument

# from patternmatcher.expression import Times, Plus, Transpose, Inverse, Identity, \
#                                       NumericConstant, \
#                                       Scalar, Vector, Matrix

# import patternmatcher.properties as properties

n = 10
m = 20
k = 30

"""
Note:
- For addition, we decided to only use a binary kernel.
"""

###############
# BLAS "0"

# scalar product

alpha = Scalar("alpha")
beta = Scalar("beta")
out = Scalar("out")
cf = lambda d: 1

scalar_product = KernelDescription(
    OperationKV(
        None,
        {None: Times(alpha, beta)}
        ),
    [], # variants
    out, # return value
    cf, # cost function
    "",
    "*$out = (*$alpha) * (*$beta);",
    "",
    [], # Argument objects
    )

# TODO scalar division?

# scalar sum

alpha = Scalar("alpha")
beta = Scalar("beta")
out = Scalar("out")
cf = lambda d: 1

scalar_sum = KernelDescription(
    OperationKV(
        None,
        {None: Plus(alpha, beta)}
    ),
    [], # variants
    out, # return value
    cf, # cost function
    "",
    "*$out = (*$alpha) + (*$beta);",
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

# TODO what about the output here?
# TODO add DOT version where on or both operands are row vectors: Different
# pattern and operands, same kernel call.

dot = KernelDescription(
    OperationKV(
        None,
        {None: Times(Transpose(x),y)}
    ),
    [], # variants
    out, # return value
    cf, # cost function
    "",
    "*$out = ${type_prefix}dot($N, $x, $incx, $y, $incy);",
    "",
    [SizeArgument("N", x, "rows"),
     StrideArgument("incx", x, "rows"),
     StrideArgument("incy", y, "rows")], # Argument objects
    )

# SCAL

x = Vector("x", (n, 1))
alpha = Scalar("alpha")
cf = lambda d: d["N"]

scal = KernelDescription(
    OperationKV(
        None,
        {None: Times(alpha, x)}
    ),
    [], # variants
    x, # return value
    cf, # cost function
    "",
    "${type_prefix}scal($N, $alpha, $x, $incx);",
    "",
    [SizeArgument("N", x, "rows"),
     StrideArgument("incx", x, "rows")], # Argument objects
    )

# AXPY

x = Vector("x", (n, 1))
y = Vector("y", (n, 1))
alpha = Scalar("alpha")
cf = lambda d: d["N"]

axpy = KernelDescription(
    OperationKV(
        None,
        {None: Plus(Times(alpha, x), y)}),
    [
        DefaultValueKV(alpha, [NumericConstant(1)])
    ],
    y, # return value
    cf, # cost function
    "",
    "${type_prefix}axpy($N, $alpha, $x, $incx, $y, $incy);",
    "",
    [SizeArgument("N", x, "rows"),
     StrideArgument("incx", x, "rows"),
     StrideArgument("incy", y, "rows")], # Argument objects
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
    OperationKV(
        None,
        {None: Plus(Times(alpha, x, Transpose(y)), A)}
    ),
    [
        DefaultValueKV(alpha, [NumericConstant(1)])
    ],
    A, # return value
    cf, # cost function
    "",
    "${type_prefix}ger($M, $N, $alpha, $x, $incx, $y, $incy, $A, $ldA);",
    "",
    [SizeArgument("M", x, "rows"),
     SizeArgument("N", y, "rows"),
     StrideArgument("incx", x, "rows"),
     StrideArgument("incy", y, "rows"),
     StrideArgument("ldA", A, "rows")], # Argument objects
    )


A = Matrix("A", (m, n))
x = Vector("x", (m, 1))
y = Vector("y", (n, 1))
alpha = Scalar("alpha")
cf = lambda d: 2*d["N"]*d["M"]

# IMPORTANT: for this (non-overwriting) case, A has to be initialized with zeros
ger_alt = KernelDescription(
    OperationKV(
        None,
        {None: Times(alpha, x, Transpose(y))}
    ),
    [
        DefaultValueKV(alpha, [NumericConstant(1)])
    ],
    A, # return value
    cf, # cost function
    "",
    "${type_prefix}ger($M, $N, $alpha, $x, $incx, $y, $incy, $A, $ldA);",
    "",
    [SizeArgument("M", x, "rows"),
     SizeArgument("N", y, "rows"),
     StrideArgument("incx", x, "rows"),
     StrideArgument("incy", y, "rows"),
     StrideArgument("ldA", A, "rows")], # Argument objects
    )


# GEMV

A = Matrix("A", (m, n))
x = Vector("x", (n, 1))
y = Vector("y", (m, 1))
alpha = Scalar("alpha")
beta = Scalar("beta")
cf = lambda d: 2*d["N"]*d["M"]

gemv = KernelDescription(
    OperationKV(
        None,
        {None: Plus(Times(alpha, Op1(A), x), Times(beta, y))}
    ),
    [
        OperatorKV("transA", {"N": Identity, "T": Transpose}, Op1),
        DefaultValueKV(alpha, [NumericConstant(1)]),
        DefaultValueKV(beta, [NumericConstant(0), NumericConstant(1)])
    ],
    y, # return value
    cf, # cost function
    "",
    "${type_prefix}gemv($transA, $M, $N, $alpha, $A, $ldA, $x, $incx, $beta, $y, $incy);",
    "",
    [SizeArgument("M", A, "rows"),
     SizeArgument("N", A, "columns"),
     StrideArgument("incx", x, "rows"),
     StrideArgument("incy", y, "rows"),
     StrideArgument("ldA", A, "rows")], # Argument objects
    [KernelType.identity, KernelType.transpose]
    )

# TRSV

A = Matrix("A", (n, n))
A.set_property(Property.SQUARE)
x = Vector("x", (n, 1))
cf = lambda d: d["N"]**2

trsv = KernelDescription(
    OperationKV(
        None,
        {None: Times(Op1(Inverse(A)), x)}
    ),
    [
        OperatorKV("transA", {"N": Identity, "T": Transpose}, Op1),
        PropertyKV("uplo", {"U": Property.UPPER_TRIANGULAR, "L": Property.LOWER_TRIANGULAR}, A)
    ],
    x, # return value
    cf, # cost function
    "",
    "${type_prefix}trsv($uplo, $transA, $diag, $N, $A, $ldA, $x, $incx);",
    "",
    [SizeArgument("N", A, "rows"),
     PropertyArgument("diag", A, Property.UNIT_DIAGONAL, ["U", "N"]),
     StrideArgument("incx", x, "rows"),
     StrideArgument("ldA", A, "rows")], # Argument objects
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

gemm = KernelDescription(
    OperationKV(
        None,
        {None: Plus(Times(alpha, Op1(A), Op2(B)), Times(beta, C))}
    ),
    [
        OperatorKV("transA", {"N": Identity, "T": Transpose}, Op1),
        OperatorKV("transB", {"N": Identity, "T": Transpose}, Op2),
        DefaultValueKV(alpha, [NumericConstant(1)]),
        DefaultValueKV(beta, [NumericConstant(0), NumericConstant(1)])
    ],
    C, # return value
    cf, # cost function
    "",
    "${type_prefix}gemm($transA, $transB, $M, $N, $K, $alpha, $A, $ldA, $B, $ldB, $beta, $C, $ldC);",
    "",
    [SizeArgument("M", Op1(A), "rows"),
     SizeArgument("N", Op2(B), "columns"),
     SizeArgument("K", Op1(A), "columns"),
     StrideArgument("ldA", A, "rows"),
     StrideArgument("ldB", B, "rows"),
     StrideArgument("ldC", C, "rows")], # Argument objects
    )

# TRSM

A = Matrix("A", (m, m))
A.set_property(Property.SQUARE)
B = Matrix("B", (m, n))
alpha = Scalar("alpha")
cf = lambda d: (d["N"]**2)*d["M"]

trsm = KernelDescription(
    OperationKV(
        "side",
        {"L": Times(alpha, Op1(Inverse(A)), B),
         "R": Times(alpha, B, Op1(Inverse(A)))}
    ),
    [
        OperatorKV("transA", {"N": Identity, "T": Transpose}, Op1),
        DefaultValueKV(alpha, [NumericConstant(1)]),
        PropertyKV("uplo", {"U": Property.UPPER_TRIANGULAR, "L": Property.LOWER_TRIANGULAR}, A)
    ],
    B, # return value
    cf, # cost function
    "",
    "${type_prefix}trsm($side, $uplo, $transA, $diag, $M, $N, $alpha, $A, $ldA, $B, $ldB);",
    "",
    [SizeArgument("M", B, "rows"),
     SizeArgument("N", B, "columns"),
     PropertyArgument("diag", A, Property.UNIT_DIAGONAL, ["U", "N"]),
     StrideArgument("ldA", A, "rows"),
     StrideArgument("ldB", B, "rows")], # Argument objects
    [KernelType.identity, KernelType.transpose]
    )

###############
# LAPACK

# LASCL (scalar * matrix) (LAPACK)
# dlascl(TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO)
# http://www.netlib.org/lapack/explore-html/d7/d43/group__aux_o_t_h_e_rauxiliary_ga7bce4c35ec5a86ee0bfdd15c476d99c8.html#ga7bce4c35ec5a86ee0bfdd15c476d99c8

A = Matrix("A", (m, n))
cfrom = Scalar("cfrom")
cto = Scalar("cto")
cf = lambda d: d["N"]*d["M"]

lascl = KernelDescription(
    OperationKV(
        None,
        {None: Times(cto, Inverse(cfrom), A)}
    ),
    [
        PropertyKV(
            "type", # TODO this should be a property argument. But that's difficult.
            {"G": Property.MATRIX, # this should be general/full/whatever. MATRIX works because there is scal for vectors.
            "U": Property.UPPER_TRIANGULAR,
            "L": Property.LOWER_TRIANGULAR}, # even more properties are possible
            A),
        DefaultValueKV(cfrom, [NumericConstant(1)]),
        DefaultValueKV(cto, [NumericConstant(1)])
    ],
    A, # return value
    cf, # cost function
    "",
    "${type_prefix}lascl($type, $kl, $ku, $cfrom, $cto, $M, $N, $A, $ldA, info);",
    "",
    [SizeArgument("M", A, "rows"),
     SizeArgument("N", A, "columns"),
    StrideArgument("ldA", A, "rows")], # Argument objects
    )

# getri (matrix inversion) (somewhat fake)
# IMPORTANT: This requires the result of an LU factorization.
# TODO how to handle this workspace thing? Solution: Each kernel has a set of
# lines of code that it depends one, some are generated, some are constant. This
# includes: compute lwork, allocate lwork memory, define constants.
# ipiv is input (expected to come from LU decomposition). To use this as a stand
# alone kernel, I need to pass the idendity 

A = Matrix("A", (n, n))
A.set_property(Property.SQUARE)
cf = lambda d: 2*d["N"]**3

getri = KernelDescription(
    OperationKV(
        None,
        {None: Inverse(A)}
    ),
    [], # variants
    A, # return value
    cf, # cost function
    """\
    $type* work_$work_id = 0;
    int lwork_$work_id = -1;
    ${type_prefix}getri($N, NULL, $ldA, NULL, $work, $lwork, info);
    lwork_$work_id = work_$work_id[0];
    work_$work_id = ($type*) malloc(lwork_$work_id);
    """,
    "${type_prefix}getri($N, $A, $ldA, $ipiv, $work, $lwork, info);",
    """free(work_$work_id);\n""",
    [SizeArgument("N", A, "columns"),
     StrideArgument("ldA", A, "rows")], # Argument objects
    [KernelType.identity, KernelType.transpose]
    )

# trtri (triangular matrix inversion)

A = Matrix("A", (n, n))
A.set_property(Property.SQUARE)
cf = lambda d: (d["N"]**3)/3

trtri = KernelDescription(
    OperationKV(
        None,
        {None: Inverse(A)}
    ),
    [
        PropertyKV("uplo", {"U": Property.UPPER_TRIANGULAR, "L": Property.LOWER_TRIANGULAR}, A)
    ],
    A, # return value
    cf, # cost function
    "",
    "${type_prefix}trtri($uplo, $diag, $N, $A, $ldA, info);",
    "",
    [SizeArgument("N", A, "rows"),
     StrideArgument("ldA", A, "rows"),
     PropertyArgument("diag", A, Property.UNIT_DIAGONAL, ["U", "N"])], # Argument objects
    [KernelType.identity, KernelType.transpose]
    )

# matrix transposition (totally fake)

A = Matrix("A", (n, n))
cf = lambda d: 0

transpose = KernelDescription(
    OperationKV(
        None,
        {None: Transpose(A)}
    ),
    [], # variants
    A, # return value
    cf, # cost function
    "",
    "transpose($N, $M, $A, $ldA);",
    "",
    [SizeArgument("M", A, "rows"),
     SizeArgument("N", A, "columns"),
     StrideArgument("ldA", A, "rows")], # Argument objects
    # [KernelType.identity, KernelType.transpose]
    )

# vector transposition (totally fake) (this is kind of tricky. Since it doesn't actually do anything, overwriting is never a problem. Solution: Just use different symbol for output?)

x = Vector("x", (n, 1))
cf = lambda d: 0

transpose_vector = KernelDescription(
    OperationKV(
        None,
        {None: Transpose(x)}
    ),
    [], # variants
    x, # return value
    cf, # cost function
    "",
    "// Transposing vector $x (no operation);",
    "",
    [], # Argument objects
    # [KernelType.identity, KernelType.transpose]
    )

# matrix sum

alpha = Scalar("alpha")
A = Matrix("A", (m, n))
B = Matrix("B", (m, n))
cf = lambda d: 2*d["N"]*d["M"]

# NOTE We are adding up columns. Thus, incx and incy are column stride
# TODO what about not doing the pointer arithmetic in the signature: i < N*ldA, i+=ldA
matrix_sum_NN = KernelDescription(
    OperationKV(
        None,
        {None: Plus(Times(alpha, A), B)} # a*x+y
    ),
    [
        DefaultValueKV(alpha, [NumericConstant(1)])
    ],
    B, # return value
    cf, # cost function
    "",
    """\
    for(int i = 0; i < $N; i++){
        ${type_prefix}axpy($MN, $alpha, $A+($rsA*i), $csA, $B+($rsB*i), $csB);
    }\
    """,
    "",
    [SizeArgument("N", A, "columns", as_value=True), # as value
     SizeArgument("M", A, "rows"), # for cost function only
     SizeArgument("MN", A, "entries"),
     StrideArgument("rsA", A, "rows", as_value=True), # as value
     StrideArgument("csA", A, "columns"),
     StrideArgument("rsB", B, "rows", as_value=True), # as value
     StrideArgument("csB", B, "columns")], # Argument objects
    )

# matrix sum with one transposition

alpha = Scalar("alpha")
A = Matrix("A", (m, n))
B = Matrix("B", (m, n))
cf = lambda d: 2*d["N"]*d["M"]

# NOTE We are adding up columns. Thus, incx and incy are column stride
# TODO what about not doing the pointer arithmetic in the signature: i < N*ldA, i+=ldA
matrix_sum_TN = KernelDescription(
    OperationKV(
    None,
    {None: Plus(Times(alpha, Transpose(A)), B)} # a*x+y
    ),
    [
        DefaultValueKV(alpha, [NumericConstant(1)])
    ],
    B, # return value
    cf, # cost function
    "",
    """\
    for(int i = 0; i < $N; i++){
        ${type_prefix}axpy($MN, $alpha, $A+($csA*i), $rsA, $B+($rsB*i), $csB);
    }\
    """,
    "",
    [SizeArgument("N", A, "columns", as_value=True), # as value
     SizeArgument("M", A, "rows"), # for cost function only
     SizeArgument("MN", A, "entries"),
     StrideArgument("rsA", A, "rows"),
     StrideArgument("csA", A, "columns", as_value=True),
     StrideArgument("rsB", B, "rows", as_value=True), # as value
     StrideArgument("csB", B, "columns")], # Argument objects
    )

# matrix sum with one transposition

alpha = Scalar("alpha")
A = Matrix("A", (m, n))
B = Matrix("B", (m, n))
cf = lambda d: 2*d["N"]*d["M"]

# NOTE We are adding up columns. Thus, incx and incy are column stride
# TODO what about not doing the pointer arithmetic in the signature: i < N*ldA, i+=ldA
matrix_sum_NT = KernelDescription(
    OperationKV(
        None,
        {None: Plus(Times(alpha, A), Transpose(B))} # a*x+y
    ),
    [
        DefaultValueKV(alpha, [NumericConstant(1)])
        # TODO remove default value of either NT or TN? With default values, they become the same (the only difference is which operand is overwritten)
    ],
    B, # return value
    cf, # cost function
    "",
    """\
    for(int i = 0; i < $N; i++){
        ${type_prefix}axpy($MN, $alpha, $A+($rsA*i), $csA, $B+($csB*i), $rsB);
    }\
    """,
    "",
    [SizeArgument("N", A, "columns", as_value=True), # as value
     SizeArgument("M", A, "rows"), # for cost function only
     SizeArgument("MN", A, "entries"),
     StrideArgument("rsA", A, "rows", as_value=True), # as value
     StrideArgument("csA", A, "columns"),
     StrideArgument("rsB", B, "rows"),
     StrideArgument("csB", B, "columns", as_value=True)], # Argument objects
    )


# matrix sum with two transpositions

alpha = Scalar("alpha")
A = Matrix("A", (m, n))
B = Matrix("B", (m, n))
cf = lambda d: 2*d["N"]*d["M"]

# NOTE We are adding up columns. Thus, incx and incy are column stride
# TODO what about not doing the pointer arithmetic in the signature: i < N*ldA, i+=ldA
matrix_sum_TT = KernelDescription(
    OperationKV(
        None,
        {None: Plus(Times(alpha, Transpose(A)), Transpose(B))} # a*x+y
    ),
    [
        DefaultValueKV(alpha, [NumericConstant(1)])
    ],
    B, # return value
    cf, # cost function
    "",
    """\
    for(int i = 0; i < $N; i++){
        ${type_prefix}axpy($MN, $alpha, $A+($rsA*i), $csA, $B+($csB*i), $rsB);
    }\
    """,
    "",
    [SizeArgument("N", A, "columns", as_value=True), # as value
     SizeArgument("M", A, "rows"), # for cost function only
     SizeArgument("MN", A, "entries"),
     StrideArgument("rsA", A, "rows"),
     StrideArgument("csA", A, "columns", as_value=True),
     StrideArgument("rsB", B, "rows"),
     StrideArgument("csB", B, "columns", as_value=True)], # Argument objects
    )