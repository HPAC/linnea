# from patternmatcher.automaton import PMA

import itertools

from . import kernels
from . import factorizations


# Kernels
other_kernel_descriptions = [
     clak.kernels.c.kernels.scalar_product,
     clak.kernels.c.kernels.scalar_sum,
     clak.kernels.c.kernels.dot,
     clak.kernels.c.kernels.ger,
     clak.kernels.c.kernels.ger_alt,
     clak.kernels.c.kernels.gemv,
     clak.kernels.c.kernels.trsv,
     clak.kernels.c.kernels.gemm,
     clak.kernels.c.kernels.trsm
     ]

others = list(itertools.chain.from_iterable( list(description.generate_kernels()) for description in other_kernel_descriptions ))

# Unary kernels
unary_kernel_descriptions = [
     clak.kernels.c.kernels.getri,
     clak.kernels.c.kernels.trtri,
     clak.kernels.c.kernels.transpose,
     clak.kernels.c.kernels.transpose_vector
     ]

unary_kernels = list(itertools.chain.from_iterable( list(description.generate_kernels()) for description in unary_kernel_descriptions ))

# Addition kernels (for AdditionGraph)
# The order is (somewhat) important. Addition comes before scaling.
addition_kernel_descriptions = [
     clak.kernels.c.kernels.axpy,
     clak.kernels.c.kernels.scal,
     clak.kernels.c.kernels.matrix_sum_NN,
     clak.kernels.c.kernels.matrix_sum_TN,
     clak.kernels.c.kernels.matrix_sum_NT,
     clak.kernels.c.kernels.matrix_sum_TT,
     clak.kernels.c.kernels.lascl,
     ]

addition_kernels = list(itertools.chain.from_iterable( list(description.generate_kernels()) for description in addition_kernel_descriptions ))

kernels = others + addition_kernels

matrix_chain_kernels = [kernel for kernel in kernels if kernel.is_matrix_chain_kernel()]

# matrix_sum = list(clak.kernels.c.kernels.matrix_sum.generate_kernels())[0]

cholesky = clak.kernels.c.factorizations.cholesky
eigendecomposition = clak.kernels.c.factorizations.eigendecomposition
ldl = clak.kernels.c.factorizations.ldl
plu = clak.kernels.c.factorizations.plu
singular_value = clak.kernels.c.factorizations.singular_value
qr_square = clak.kernels.c.factorizations.qr_square
qr_column = clak.kernels.c.factorizations.qr_column

factorizations = (cholesky, eigendecomposition, ldl, plu, singular_value, qr_square, qr_column)
# factorizations = (cholesky, eigendecomposition, ldl, plu, QR())
factorizations_by_type = ((cholesky, ldl, plu), (qr_square, qr_column), (eigendecomposition, singular_value))


MC_automaton = PMA()
MC_automaton.add_patterns([kernel.pattern for kernel in matrix_chain_kernels])

