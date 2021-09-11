
import itertools

from . import kernels
from . import factorizations

import matchpy

def to_kernels(descriptions):
    return list(itertools.chain.from_iterable( list(description.generate_kernels()) for description in descriptions))

scalar_product = list(kernels.scalar_product.generate_kernels())
scalar_division = list(kernels.scalar_division.generate_kernels())
scalar_inversion = list(kernels.scalar_inversion.generate_kernels())
scalar_sum = list(kernels.scalar_sum.generate_kernels())
dot = list(kernels.dot.generate_kernels())
ger = list(kernels.ger.generate_kernels())
ger_alt = list(kernels.ger_alt.generate_kernels())
gemv = list(kernels.gemv.generate_kernels())
trsv = list(kernels.trsv.generate_kernels())
gemm = list(kernels.gemm.generate_kernels())
trsm = list(kernels.trsm.generate_kernels())
getri = list(kernels.getri.generate_kernels())
trtri = list(kernels.trtri.generate_kernels())
syinv = list(kernels.syinv.generate_kernels())
diaginv = list(kernels.diaginv.generate_kernels())
transpose = list(kernels.transpose.generate_kernels())
transpose_vector = list(kernels.transpose_vector.generate_kernels())
axpy = list(kernels.axpy.generate_kernels())
scal = list(kernels.scal.generate_kernels())
matrix_sum = list(kernels.matrix_sum.generate_kernels())
matrix_sum_transpose = list(kernels.matrix_sum_transpose.generate_kernels())
lascl = list(kernels.lascl.generate_kernels())
syrk = list(kernels.syrk.generate_kernels())
syr2k = list(kernels.syr2k.generate_kernels())
symv = list(kernels.symv.generate_kernels())
trmv = list(kernels.trmv.generate_kernels())
syr = list(kernels.syr.generate_kernels())
syr_alt = list(kernels.syr_alt.generate_kernels())
trmm = list(kernels.trmm.generate_kernels())
symm = list(kernels.symm.generate_kernels())
diagdiagmul = list(kernels.diagdiagmul.generate_kernels())
diagdiagsolve = list(kernels.diagdiagsolve.generate_kernels())
diagsmr = list(kernels.diagsmr.generate_kernels())
diagsml = list(kernels.diagsml.generate_kernels())
diagsv = list(kernels.diagsv.generate_kernels())
diagmmr = list(kernels.diagmmr.generate_kernels())
diagmml = list(kernels.diagmml.generate_kernels())
diagmv = list(kernels.diagmv.generate_kernels())
diagscal = list(kernels.diagscal.generate_kernels())
diagdiagadd = list(kernels.diagdiagadd.generate_kernels())
diagfulladd = list(kernels.diagfulladd.generate_kernels())
invscal = list(kernels.invscal.generate_kernels())
invlascl = list(kernels.invlascl.generate_kernels())
# direct solvers (so far for matrix chain paper only. Problem: They overwrite A and B, which can not be expressed at the moment)
posv = list(kernels.posv.generate_kernels())
posvr = list(kernels.posvr.generate_kernels())
sysv = list(kernels.sysv.generate_kernels())
sysvr = list(kernels.sysvr.generate_kernels())
gesv = list(kernels.gesv.generate_kernels())
gesvr = list(kernels.gesvr.generate_kernels())
gesvrt = list(kernels.gesvrt.generate_kernels())
posv_vec = list(kernels.posv_vec.generate_kernels())
sysv_vec = list(kernels.sysv_vec.generate_kernels())
gesv_vec = list(kernels.gesv_vec.generate_kernels())
# operations with permutation matrices
pmm = list(kernels.pmm.generate_kernels())
ptmm = list(kernels.ptmm.generate_kernels())
mpm = list(kernels.mpm.generate_kernels())
mptm = list(kernels.mptm.generate_kernels())
pvm = list(kernels.pvm.generate_kernels())
ptvm = list(kernels.ptvm.generate_kernels())
ppm = list(kernels.ppm.generate_kernels())
ptpm = list(kernels.ptpm.generate_kernels())
pptm = list(kernels.pptm.generate_kernels())
transpose_perm = list(kernels.transpose_perm.generate_kernels())

# Currently, this list is not used (except for obtaining matrix chain kernels.)
kernels = list(itertools.chain(
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
                matrix_sum_transpose,
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
                ## Solvers, used for CGO 2018.
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
                ## end
                pmm,
                ptmm,
                mpm,
                mptm,
                pvm,
                ptvm,
                ppm,
                ptpm,
                pptm,
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
                matrix_sum_transpose,
                diagdiagadd,
                diagfulladd,
            ))

scaling_kernels = list(itertools.chain(
                scal,
                lascl,
                diagscal,
                # It's debatable whether scalar_product belongs here. It's here
                # to make decompose_sum work for something like alpha+2*beta.
                scalar_product,  
            ))


matrix_chain_kernels = [kernel for kernel in kernels if kernel.is_matrix_chain_kernel()]

cholesky = factorizations.cholesky
eigendecomposition = factorizations.eigendecomposition
plu = factorizations.plu
singular_value_sq = factorizations.singular_value_sq
singular_value_cp = factorizations.singular_value_cp
singular_value_rp = factorizations.singular_value_rp
qr_square = factorizations.qr_square
qr_column = factorizations.qr_column

factorizations = (cholesky, eigendecomposition, plu, singular_value_sq, singular_value_cp, singular_value_rp, qr_square, qr_column)
factorizations_by_type = ((cholesky, plu), (qr_square, qr_column), (eigendecomposition, singular_value_sq, singular_value_cp, singular_value_rp))

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

kernel_MA = matchpy.ManyToOneMatcher()
for kernel in kernels:
    kernel_MA.add(kernel.pattern_with_context, label=kernel)

def clear():
    global addition_kernel_MA, kernel_MA
    addition_kernel_MA.clear()
    kernel_MA.clear()
