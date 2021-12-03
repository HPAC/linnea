using LinearAlgebra.BLAS
using LinearAlgebra

"""
    algorithm4(ml0::Array{Float64,2}, ml1::Array{Float64,1})

Compute
b = ((X^T X)^-1 X^T y).

Requires at least Julia v1.0.

# Arguments
- `ml0::Array{Float64,2}`: Matrix X of size 1500 x 1000 with property FullRank.
- `ml1::Array{Float64,1}`: Vector y of size 1500.
"""                    
function algorithm4(ml0::Array{Float64,2}, ml1::Array{Float64,1})
    # cost: 1.84e+09 FLOPs
    # X: ml0, full, y: ml1, full
    ml2 = Array{Float64}(undef, 1000)
    # tmp8 = (X^T y)
    gemv!('T', 1.0, ml0, ml1, 0.0, ml2)

    # X: ml0, full, tmp8: ml2, full
    ml3 = Array{Float64}(undef, 1000, 1000)
    # tmp1 = (X^T X)
    syrk!('L', 'T', 1.0, ml0, 0.0, ml3)

    # tmp8: ml2, full, tmp1: ml3, symmetric_lower_triangular
    # (L2 L2^T) = tmp1
    LAPACK.potrf!('L', ml3)

    # tmp8: ml2, full, L2: ml3, lower_triangular
    # tmp10 = (L2^-1 tmp8)
    trsv!('L', 'N', 'N', ml3, ml2)

    # L2: ml3, lower_triangular, tmp10: ml2, full
    # tmp11 = (L2^-T tmp10)
    trsv!('L', 'T', 'N', ml3, ml2)

    # tmp11: ml2, full
    # b = tmp11
    return (ml2)
end