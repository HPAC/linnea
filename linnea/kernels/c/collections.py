# from patternmatcher.automaton import PMA

import itertools

from . import reductions
from . import factorizations


# Normal reduction kernels
other_kernel_descriptions = [
     clak.kernels.c.reductions.scalar_product,
     clak.kernels.c.reductions.scalar_sum,
     clak.kernels.c.reductions.dot,
     clak.kernels.c.reductions.ger,
     clak.kernels.c.reductions.ger_alt,
     clak.kernels.c.reductions.gemv,
     clak.kernels.c.reductions.trsv,
     clak.kernels.c.reductions.gemm,
     clak.kernels.c.reductions.trsm
     ]

others = list(itertools.chain.from_iterable( list(description.generate_kernels()) for description in other_kernel_descriptions ))

# Unary reduction kernels
unary_kernel_descriptions = [
     clak.kernels.c.reductions.getri,
     clak.kernels.c.reductions.trtri,
     clak.kernels.c.reductions.transpose,
     clak.kernels.c.reductions.transpose_vector
     ]

unary_kernels = list(itertools.chain.from_iterable( list(description.generate_kernels()) for description in unary_kernel_descriptions ))

# Addition kernels (for AdditionGraph)
# The order is (somewhat) important. Addition comes before scaling.
addition_kernel_descriptions = [
     clak.kernels.c.reductions.axpy,
     clak.kernels.c.reductions.scal,
     clak.kernels.c.reductions.matrix_sum_NN,
     clak.kernels.c.reductions.matrix_sum_TN,
     clak.kernels.c.reductions.matrix_sum_NT,
     clak.kernels.c.reductions.matrix_sum_TT,
     clak.kernels.c.reductions.lascl,
     ]

addition_kernels = list(itertools.chain.from_iterable( list(description.generate_kernels()) for description in addition_kernel_descriptions ))

reductions = others + addition_kernels

matrix_chain_kernels = [kernel for kernel in reductions if kernel.is_matrix_chain_kernel()]

# matrix_sum = list(clak.kernels.c.reductions.matrix_sum.generate_kernels())[0]

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

