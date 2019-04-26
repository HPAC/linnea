
import itertools

from . import reductions
from . import factorizations

import matchpy

def to_kernels(descriptions):
    return list(itertools.chain.from_iterable( list(description.generate_kernels()) for description in descriptions))

scalar_product = list(reductions.scalar_product.generate_kernels())
scalar_division = list(reductions.scalar_division.generate_kernels())
scalar_inversion = list(reductions.scalar_inversion.generate_kernels())
scalar_sum = list(reductions.scalar_sum.generate_kernels())
dot = list(reductions.dot.generate_kernels())
ger = list(reductions.ger.generate_kernels())
ger_alt = list(reductions.ger_alt.generate_kernels())
gemv = list(reductions.gemv.generate_kernels())
trsv = list(reductions.trsv.generate_kernels())
gemm = list(reductions.gemm.generate_kernels())
trsm = list(reductions.trsm.generate_kernels())
getri = list(reductions.getri.generate_kernels())
trtri = list(reductions.trtri.generate_kernels())
syinv = list(reductions.syinv.generate_kernels())
diaginv = list(reductions.diaginv.generate_kernels())
transpose = list(reductions.transpose.generate_kernels())
transpose_vector = list(reductions.transpose_vector.generate_kernels())
axpy = list(reductions.axpy.generate_kernels())
scal = list(reductions.scal.generate_kernels())
matrix_sum = list(reductions.matrix_sum.generate_kernels())
lascl = list(reductions.lascl.generate_kernels())
syrk = list(reductions.syrk.generate_kernels())
syr2k = list(reductions.syr2k.generate_kernels())
symv = list(reductions.symv.generate_kernels())
trmv = list(reductions.trmv.generate_kernels())
syr = list(reductions.syr.generate_kernels())
syr_alt = list(reductions.syr_alt.generate_kernels())
trmm = list(reductions.trmm.generate_kernels())
symm = list(reductions.symm.generate_kernels())
diagdiagmul = list(reductions.diagdiagmul.generate_kernels())
diagdiagsolve = list(reductions.diagdiagsolve.generate_kernels())
diagsmr = list(reductions.diagsmr.generate_kernels())
diagsml = list(reductions.diagsml.generate_kernels())
diagsv = list(reductions.diagsv.generate_kernels())
diagmmr = list(reductions.diagmmr.generate_kernels())
diagmml = list(reductions.diagmml.generate_kernels())
diagmv = list(reductions.diagmv.generate_kernels())
diagscal = list(reductions.diagscal.generate_kernels())
diagdiagadd = list(reductions.diagdiagadd.generate_kernels())
diagfulladd = list(reductions.diagfulladd.generate_kernels())
invscal = list(reductions.invscal.generate_kernels())
invlascl = list(reductions.invlascl.generate_kernels())
# direct solvers (so far for matrix chain paper only. Problem: They overwrite A and B, which can not be expressed at the moment)
posv = list(reductions.posv.generate_kernels())
posvr = list(reductions.posvr.generate_kernels())
sysv = list(reductions.sysv.generate_kernels())
sysvr = list(reductions.sysvr.generate_kernels())
gesv = list(reductions.gesv.generate_kernels())
gesvr = list(reductions.gesvr.generate_kernels())
gesvrt = list(reductions.gesvrt.generate_kernels())
posv_vec = list(reductions.posv_vec.generate_kernels())
sysv_vec = list(reductions.sysv_vec.generate_kernels())
gesv_vec = list(reductions.gesv_vec.generate_kernels())
# operations with permutation matrices
pmm = list(reductions.pmm.generate_kernels())
ptmm = list(reductions.ptmm.generate_kernels())
mpm = list(reductions.mpm.generate_kernels())
mptm = list(reductions.mptm.generate_kernels())
pvm = list(reductions.pvm.generate_kernels())
ptvm = list(reductions.ptvm.generate_kernels())
ppm = list(reductions.ppm.generate_kernels())
ptpm = list(reductions.ptpm.generate_kernels())
transpose_perm = list(reductions.transpose_perm.generate_kernels())

# Currently, this list is not used (except for obtaining matrix chain kernels.)
reductions = list(itertools.chain(
                scalar_product,
                scalar_division,
                scalar_sum,
                dot,
                ger,
                ger_alt,
                gemv,
                trsv,
                gemm,
                trsm,
                axpy,
                scal,
                matrix_sum,
                lascl,
                symm,
                syrk,
                syr2k,
                symv,
                trmv,
                syr,
                syr_alt,
                trmm,
                diagscal,
                diagdiagmul,
                diagdiagsolve,
                diagsmr,
                diagsml,
                diagsv,
                diagmmr,
                diagmml,
                diagmv,
                diagdiagadd,
                diagfulladd,
                invscal,
                invlascl,
                # posv,
                # posvr,
                # sysv,
                # sysvr,
                # gesv,
                # gesvr,
                # gesvrt,
                # posv_vec,
                # sysv_vec,
                # gesv_vec,
                pmm,
                ptmm,
                mpm,
                mptm,
                pvm,
                ptvm,
                ppm,
                ptpm,
            ))

unused = list(itertools.chain(
            ))

# inversion_kernels = list(itertools.chain(
#                 getri,
#                 trtri,
#                 diaginv,
#             ))

transposition_kernels = list(itertools.chain(
                transpose,
                transpose_vector,
                transpose_perm
            ))

unary_kernels = list(itertools.chain(
                # getri,
                trtri,
                diaginv,
                transpose,
                transpose_vector,
                transpose_perm,
                scalar_inversion,
                # syinv # this is horribly slow
            ))

addition_kernels = list(itertools.chain(
                scalar_sum,
                axpy,
                matrix_sum,
                diagdiagadd,
                diagfulladd,
            ))

scaling_kernels = list(itertools.chain(
                scal,
                lascl,
                diagscal,
            ))


matrix_chain_kernels = [kernel for kernel in reductions if kernel.is_matrix_chain_kernel()]

cholesky = factorizations.cholesky
eigendecomposition = factorizations.eigendecomposition
# ldl = factorizations.ldl
plu = factorizations.plu
singular_value = factorizations.singular_value
qr_square = factorizations.qr_square
qr_column = factorizations.qr_column

# factorizations = (cholesky, eigendecomposition, ldl, plu, singular_value, qr_square, qr_column)
factorizations = (cholesky, eigendecomposition, plu, singular_value, qr_square, qr_column)
# factorizations_by_type = ((cholesky, ldl, plu), (qr_square, qr_column), (eigendecomposition, singular_value))
factorizations_by_type = ((cholesky, plu), (qr_square, qr_column), (eigendecomposition, singular_value))

matrix_chain_DN = matchpy.DiscriminationNet()
for kernel in matrix_chain_kernels:
    matrix_chain_DN.add(kernel.pattern, final_label=kernel)

unary_kernel_DN = matchpy.DiscriminationNet()
for kernel in unary_kernels:
    unary_kernel_DN.add(kernel.pattern, final_label=kernel)

transposition_kernel_DN = matchpy.DiscriminationNet()
for kernel in transposition_kernels:
    transposition_kernel_DN.add(kernel.pattern, final_label=kernel)

scaling_kernel_DN = matchpy.DiscriminationNet()
for kernel in scaling_kernels:
    scaling_kernel_DN.add(kernel.pattern, final_label=kernel)

addition_kernel_MA = matchpy.ManyToOneMatcher()
for kernel in addition_kernels:
    addition_kernel_MA.add(kernel.pattern_with_context, label=kernel)

reduction_MA = matchpy.ManyToOneMatcher()
for kernel in reductions:
    reduction_MA.add(kernel.pattern_with_context, label=kernel)
    
