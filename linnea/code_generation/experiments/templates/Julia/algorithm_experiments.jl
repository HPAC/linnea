using LinearAlgebra.BLAS
using LinearAlgebra

function {}({})
    start::Float64 = 0.0
    finish::Float64 = 0.0
    Benchmarker.cachescrub()
    Benchmarker.gcscrub()
    start = time_ns()

{}

    finish = time_ns()
    return (tuple({}), (finish-start)*1e-9)
end