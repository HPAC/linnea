from ..algebra.expression import Symbol, Scalar, Vector, Matrix, ConstantScalar, \
                                 Equal, Plus, Times, Transpose, Inverse, \
                                 InverseTranspose, InverseConjugate, \
                                 InverseConjugateTranspose, \
                                 ConjugateTranspose, Index, IdentityMatrix

from ..algebra.properties import Property as properties

from ..algebra.equations import Equations

from .. import derivation



class Example01():
    def __init__(self):

        # least squares

        n = 1500
        m = 1000

        X = Matrix("X", (n, m))
        X.set_property(properties.INPUT)
        X.set_property(properties.FULL_RANK)

        y = Vector("y", (n, 1))
        y.set_property(properties.INPUT)

        b = Vector("b", (m, 1))
        b.set_property(properties.OUTPUT)

        # b = (X^T X)^-1 X^T y
        self.eqns = Equations(Equal(b, Times(Inverse(Times(Transpose(X), X)), Transpose(X), y)))


class Example02():
    def __init__(self):

        # generalized least squares

        n = 1500
        m = 1000

        S = Matrix("S", (n, n))
        S.set_property(properties.SPD)
        S.set_property(properties.INPUT)

        X = Matrix("X", (n, m))
        X.set_property(properties.FULL_RANK)
        X.set_property(properties.INPUT)

        z = Vector("z", (m, 1))
        z.set_property(properties.OUTPUT)

        y = Vector("y", (n, 1))
        y.set_property(properties.INPUT)

        # clak.special_properties.add_expression(Times(Transpose(X), Inverse(S), X ), set([properties.SPD))

        self.eqns = Equations(Equal(z, Times(Inverse(Times(Transpose(X), Inverse(S), X ) ), Transpose(X), Inverse(S), y)))


class Example03():
    def __init__(self, ):

        # optimization problem 1

        n = 1500
        m = 1000

        A = Matrix("A", (m, n))
        A.set_property(properties.FULL_RANK)
        A.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        # W is positive
        W.set_property(properties.FULL_RANK)
        W.set_property(properties.DIAGONAL)
        W.set_property(properties.SPD)
        W.set_property(properties.INPUT)

        b = Vector("b", (m, 1))
        b.set_property(properties.INPUT)

        c = Vector("c", (n, 1))
        c.set_property(properties.INPUT)

        x = Vector("x", (n, 1))
        x.set_property(properties.OUTPUT)

        minus1 = ConstantScalar(-1.0)


        self.eqns = Equations(
                            Equal(
                                x,
                                Times(
                                    W,
                                    Plus(
                                        Times(Transpose(A), Inverse(Times(A, W, Transpose(A))), b),
                                        Times(minus1, c)
                                        )
                                    )
                                )
                            )


class Example04():
    def __init__(self, ):

        # optimization problem 2

        n = 1500
        m = 1000

        A = Matrix("A", (m, n))
        A.set_property(properties.FULL_RANK)
        A.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        # W is positive
        W.set_property(properties.FULL_RANK)
        W.set_property(properties.DIAGONAL)
        W.set_property(properties.SPD)
        W.set_property(properties.INPUT)

        b = Vector("b", (m, 1))
        b.set_property(properties.INPUT)

        c = Vector("c", (n, 1))
        c.set_property(properties.INPUT)

        x = Vector("x", (n, 1))
        x.set_property(properties.INPUT)

        xf = Vector("xf", (n, 1))
        xf.set_property(properties.OUTPUT)

        xo = Vector("xo", (n, 1))
        xo.set_property(properties.OUTPUT)

        minus1 = ConstantScalar(-1.0)


        
        self.eqns = Equations(
                        Equal(xf, 
                            Times(
                                W,
                                Transpose(A),
                                Inverse(Times(A, W, Transpose(A))),
                                Plus(b, Times(minus1, A, x))
                            )
                        ),
                        Equal(xo,
                            Times(W,
                                Plus(
                                    Times(
                                        Transpose(A),
                                        Inverse(Times(A, W, Transpose(A))),
                                        A, x),
                                    Times(minus1, c)
                                )
                            )
                        )
                    )


class Example05():
    def __init__(self, N = 1000):

        # signal processing

        A = Matrix("A", (N, N), properties = [properties.FULL_RANK, properties.INPUT])
        # A - tridiagonal, full rank, something even more specific
        B = Matrix("B", (N, N), properties = [properties.FULL_RANK, properties.INPUT])
        # B - tridiagonal, full rank, something even more specific
        R = Matrix("R", (N - 1, N), properties = [properties.FULL_RANK, properties.UPPER_TRIANGULAR, properties.INPUT])
        # R - upper bidiagonal
        L = Matrix("L", (N - 1, N - 1), properties = [properties.FULL_RANK, properties.DIAGONAL, properties.INPUT])

        y = Vector("y", (N, 1), properties = [properties.INPUT])
        x = Vector("x", (N, 1), properties = [properties.OUTPUT])

        self.eqns = Equations(
                            Equal(
                                x,
                                Times(
                                    Inverse(
                                        Plus(
                                            Times(
                                                InverseTranspose(A),
                                                Transpose(B),
                                                B,
                                                Inverse(A)
                                            ),
                                            Times(
                                                Transpose(R),
                                                L,
                                                R
                                            )
                                        )
                                    ),
                                    InverseTranspose(A),
                                    Transpose(B),
                                    B,
                                    Inverse(A),
                                    y
                                )
                            )
                    )


class Example06():
    def __init__(self):

        # inversion of lower triangular matrix

        n = 1000
        m = 1000
        k = 1000

        L00 = Matrix("L00", (n, n), properties = [properties.FULL_RANK, properties.LOWER_TRIANGULAR, properties.INPUT])
        L11 = Matrix("L11", (m, m), properties = [properties.FULL_RANK, properties.LOWER_TRIANGULAR, properties.INPUT])
        L22 = Matrix("L22", (k, k), properties = [properties.FULL_RANK, properties.LOWER_TRIANGULAR, properties.INPUT])
        L21 = Matrix("L21", (k, m), properties = [properties.FULL_RANK, properties.INPUT])
        L10 = Matrix("L10", (m, n), properties = [properties.FULL_RANK, properties.INPUT])
        L20 = Matrix("L20", (k, n), properties = [properties.FULL_RANK, properties.INPUT])

        X21 = Matrix("X21", (k, m), properties = [properties.OUTPUT])
        X11 = Matrix("X11", (m, m), properties = [properties.OUTPUT])
        X10 = Matrix("X10", (m, n), properties = [properties.OUTPUT])
        
        X20 = Matrix("X20", (k, n), properties = [properties.OUTPUT])
        minus1 = ConstantScalar(-1.0)

        self.eqns = Equations(
                            Equal(
                                X10,
                                Times(
                                    L10,
                                    Inverse(
                                        L00
                                    )
                                )
                            ),
                            Equal(
                                X20,
                                Plus(
                                    L20,
                                    Times(
                                        Inverse(L22),
                                        L21,
                                        Inverse(L11),
                                        L10
                                    )
                                )
                            ),
                            Equal(
                                X11,
                                Inverse(L11)
                            ),
                            Equal(
                                X21,
                                Times(
                                    minus1,
                                    Inverse(L22),
                                    L21
                                )
                            )
                    )


class Example07():
    def __init__(self):

        # local assimilation for parallel ensemble Kalman filter based on modified Cholesky decomposition 1

        N = 1000
        msd = 1000
        nsd = 1000

        minus1 = ConstantScalar(-1.0)

        B = Matrix("B", (N, N), properties = [properties.FULL_RANK, properties.SYMMETRIC, properties.INPUT])
        # B is a covariance matrix (symmetric positive semi-definite)
        H = Matrix("H", (msd, N), properties = [properties.FULL_RANK, properties.INPUT])
        R = Matrix("R", (msd, msd), properties = [properties.FULL_RANK, properties.SYMMETRIC, properties.INPUT])
        # R is a covariance matrix (symmetric positive semi-definite)
        Y = Matrix("Y", (msd, N), properties = [properties.FULL_RANK, properties.INPUT])
        Xb = Matrix("Xb", (nsd, N), properties = [properties.FULL_RANK, properties.INPUT])
        Xa = Matrix("Xa", (nsd, N), properties = [properties.FULL_RANK, properties.OUTPUT])

        self.eqns = Equations(
                            Equal(
                                Xa,
                                Plus(
                                    Xb,
                                    Times(
                                        Inverse(Plus(Inverse(B), Times(Transpose(H), Inverse(R), H))),
                                        Plus(Y, Times(minus1, H, Xb))
                                        )
                                    )
                                )
                            )


class Example08():
    def __init__(self):

        # local assimilation for parallel ensemble Kalman filter based on modified Cholesky decomposition 2

        N = 1000
        msd = 1000
        nsd = 1000

        minus1 = ConstantScalar(-1.0)

        B = Matrix("B", (N, N), properties = [properties.FULL_RANK, properties.SYMMETRIC, properties.INPUT])
        # B is a covariance matrix (symmetric positive semi-definite)
        H = Matrix("H", (msd, N), properties = [properties.FULL_RANK, properties.INPUT])
        R = Matrix("R", (msd, msd), properties = [properties.FULL_RANK, properties.SYMMETRIC, properties.INPUT])
        # R is a covariance matrix (symmetric positive semi-definite)
        Y = Matrix("Y", (msd, N), properties = [properties.FULL_RANK, properties.INPUT])
        Xb = Matrix("Xb", (nsd, N), properties = [properties.FULL_RANK, properties.INPUT])
        dX = Matrix("dX", (nsd, N), properties = [properties.FULL_RANK, properties.OUTPUT])

        self.eqns = Equations(
                        Equal(dX,
                            Times(
                                Inverse(Plus(Inverse(B), Times(Transpose(H), Inverse(R), H))),
                                Transpose(H),
                                Inverse(R),
                                Plus(Y, Times(minus1, H, Xb))
                                )
                            )
                        )

class Example09():
    def __init__(self):

        # local assimilation for parallel ensemble Kalman filter based on modified Cholesky decomposition 3

        N = 1000
        msd = 1000
        nsd = 1000

        minus1 = ConstantScalar(-1.0)

        B = Matrix("B", (N, N), properties = [properties.FULL_RANK, properties.SYMMETRIC, properties.INPUT])
        # B is a covariance matrix (symmetric positive semi-definite)
        H = Matrix("H", (msd, N), properties = [properties.FULL_RANK, properties.INPUT])
        R = Matrix("R", (msd, msd), properties = [properties.FULL_RANK, properties.SYMMETRIC, properties.INPUT])
        # R is a covariance matrix (symmetric positive semi-definite)
        Y = Matrix("Y", (msd, N), properties = [properties.FULL_RANK, properties.INPUT])
        X = Matrix("X", (nsd, N), properties = [properties.FULL_RANK, properties.INPUT])
        Xb = Matrix("Xb", (nsd, N), properties = [properties.FULL_RANK, properties.INPUT])
        dX = Matrix("dX", (nsd, N), properties = [properties.FULL_RANK, properties.OUTPUT])

        self.eqns = Equations(
                        Equal(dX,
                            Times(
                                X,
                                Transpose(Times(H, X)),
                                Inverse(Plus(R, Times(H, X, Transpose(Times(H, X))))),
                                Plus(Y, Times(minus1, H, Xb))
                                )
                            )
                        )


class Example10():
    def __init__(self):

        # image restoration 1

        n = 1500
        m = 1000

        minus1 = ConstantScalar(-1.0)
        lambda_ = Scalar("lambda") # positive
        sigma_ = Scalar("sigma_sq") # positive

        H = Matrix("H", (m, n), properties = [properties.INPUT, properties.FULL_RANK])
        I = IdentityMatrix(n, n)

        v_k = Vector("v_k", (n, 1), properties = [properties.INPUT])
        u_k = Vector("u_k", (n, 1), properties = [properties.INPUT])
        y = Vector("y", (m, 1), properties = [properties.INPUT])
        x = Vector("x", (n, 1), properties = [properties.OUTPUT])


        self.init = lambda: derivation.special_properties.add_expression(Plus(
                                        Times(
                                            Transpose(H),
                                            H
                                        ),
                                        Times(
                                            lambda_,
                                            sigma_,
                                            I
                                        )
                                    ), {properties.SPD})
        self.init()

        self.eqns = Equations(
                            Equal(
                                x,
                                Times(
                                    # (H^t * H + lambda * sigma^2 * I_n)^-1
                                    Inverse( Plus(
                                        Times(
                                            Transpose(H),
                                            H
                                        ),
                                        Times(
                                            lambda_,
                                            sigma_,
                                            I
                                        )
                                    )),
                                    # (H^T * y + lambda * sigma^2 * (v - u))
                                    Plus(
                                        Times(
                                            Transpose(H),
                                            y
                                        ),
                                        Times(
                                            lambda_,
                                            sigma_,
                                            Plus(
                                                v_k,
                                                Times(
                                                    minus1,
                                                    u_k
                                                )
                                            )
                                        )
                                    )
                                )
                            )
                        )


class Example11():
    def __init__(self):

        # image restoration 2

        n = 1500
        m = 1000

        minus1 = ConstantScalar(-1.0)
        lambda_ = Scalar("lambda") # positive
        sigma_ = Scalar("sigma_sq") # positive

        H = Matrix("H", (m, n), properties = [properties.INPUT, properties.FULL_RANK])
        H_dag_I = Matrix("H_dag_I", (n, m), properties = [properties.FULL_RANK, properties.INPUT])
        H_dag_O = Matrix("H_dag_O", (n, m), properties = [properties.OUTPUT])
        I = IdentityMatrix(n, n)

        y_k = Vector("y_k", (n, 1), properties = [properties.OUTPUT])
        y = Vector("y", (m, 1), properties = [properties.INPUT])
        x = Vector("x", (n, 1), properties = [properties.INPUT])

        h_dag = Times(
                    Transpose(H),
                    Inverse(
                        Times(
                            H,
                            Transpose(H)
                        )
                    )
                )

        self.eqns = Equations(
                        Equal(
                            H_dag_O,
                            h_dag
                        ),
                        Equal(
                            y_k,
                            Plus(
                                Times(
                                    H_dag_I,
                                    y
                                ),
                                Times(
                                    Plus(
                                        I,
                                        Times(
                                            minus1,
                                            H_dag_I,
                                            H
                                        )
                                    ),
                                    x
                                )
                            )
                        )
                    )


class Example12():
    def __init__(self):

        # image restoration 3

        n = 1500
        m = 1000

        minus1 = ConstantScalar(-1.0)
        lambda_ = Scalar("lambda") # positive
        sigma_ = Scalar("sigma_sq") # positive

        H = Matrix("H", (m, n), properties = [properties.INPUT, properties.FULL_RANK])
        H_dag_I = Matrix("H_dag_I", (n, m), properties = [properties.FULL_RANK, properties.INPUT])
        H_dag_O = Matrix("H_dag_O", (n, m), properties = [properties.OUTPUT])
        I = IdentityMatrix(n, n)

        y_k = Vector("y_k", (n, 1), properties = [properties.OUTPUT])
        y = Vector("y", (m, 1), properties = [properties.INPUT])
        x = Vector("x", (n, 1), properties = [properties.INPUT])

        h_dag = Times(
                    Transpose(H),
                    Inverse(
                        Times(
                            H,
                            Transpose(H)
                        )
                    )
                )

        self.eqns = Equations(
                        Equal(
                            y_k,
                            Plus(
                                Times(
                                    h_dag,
                                    y
                                ),
                                Times(
                                    Plus(
                                        I,
                                        Times(
                                            minus1,
                                            h_dag,
                                            H
                                        )
                                    ),
                                    x
                                )
                            )
                        )
                    )


class Example13():
    def __init__(self):

        # randomized matrix inversion 1

        # q << n
        n = 5000
        q = 500

        W = Matrix("W", (n, n))
        W.set_property(properties.SPD)
        W.set_property(properties.INPUT)

        S = Matrix("S", (n, q))
        S.set_property(properties.FULL_RANK)
        S.set_property(properties.INPUT)

        A = Matrix("A", (n, n))
        A.set_property(properties.FULL_RANK)
        A.set_property(properties.INPUT)

        Xin = Matrix("Xin", (n, n))
        Xin.set_property(properties.FULL_RANK)
        Xin.set_property(properties.INPUT)

        Xout = Matrix("Xout", (n, n))
        Xout.set_property(properties.FULL_RANK)
        Xout.set_property(properties.OUTPUT)

        I = IdentityMatrix(n, n)
        minus1 = ConstantScalar(-1.0)

        self.eqns = Equations(
                        Equal(Xout,
                            Plus(Xin,
                                Times(
                                    W, Transpose(A), S,
                                    Inverse(Times(Transpose(S), A, W, Transpose(A), S)),
                                    Transpose(S),
                                    Plus(I, Times(minus1, A, Xin))
                                    )
                                )
                            )
                        )

class Example14():
    def __init__(self):

        # randomized matrix inversion 2

        # q << n
        n = 5000
        q = 500

        W = Matrix("W", (n, n))
        W.set_property(properties.SPD)
        W.set_property(properties.INPUT)

        S = Matrix("S", (n, q))
        S.set_property(properties.FULL_RANK)
        S.set_property(properties.INPUT)

        A = Matrix("A", (n, n))
        A.set_property(properties.FULL_RANK)
        A.set_property(properties.INPUT)

        Xin = Matrix("Xin", (n, n))
        Xin.set_property(properties.FULL_RANK)
        Xin.set_property(properties.INPUT)

        Xout = Matrix("Xout", (n, n))
        Xout.set_property(properties.FULL_RANK)
        Xout.set_property(properties.OUTPUT)

        I = IdentityMatrix(n, n)
        minus1 = ConstantScalar(-1.0)

        self.eqns = Equations(
                        Equal(Xout,
                            Plus(Xin,
                                Times(
                                    Plus(I, Times(minus1, Xin, Transpose(A))),
                                    S,
                                    Inverse(Times(Transpose(S), Transpose(A), W, A, S)),
                                    Transpose(S), Transpose(A), W
                                    )
                                )
                            )
                        )


class Example15():
    def __init__(self):

        # randomized matrix inversion 3

        # q << n
        n = 4000
        q = 400

        W = Matrix("W", (n, n))
        W.set_property(properties.SPD)
        W.set_property(properties.INPUT)

        S = Matrix("S", (n, q))
        S.set_property(properties.FULL_RANK)
        S.set_property(properties.INPUT)

        A = Matrix("A", (n, n))
        A.set_property(properties.SYMMETRIC)
        A.set_property(properties.FULL_RANK)
        A.set_property(properties.INPUT)

        Xin = Matrix("Xin", (n, n))
        Xin.set_property(properties.FULL_RANK)
        Xin.set_property(properties.INPUT)

        Xout = Matrix("Xout", (n, n))
        Xout.set_property(properties.FULL_RANK)
        Xout.set_property(properties.OUTPUT)

        I = IdentityMatrix(n, n)
        minus1 = ConstantScalar(-1.0)

        Lambda = Times(S, Inverse(Times(Transpose(S), A, W, A, S)), Transpose(S))
        Omega = Times(Lambda, A, W)
        M_k = Plus(Times(Xin, A), Times(minus1, I)) 

        self.eqns = Equations(
                        Equal(Xout,
                            Plus(Xin,
                                Times(minus1, M_k, Omega),
                                Times(minus1, Transpose(Times(M_k, Omega))),
                                Times(
                                    Transpose(Omega),
                                    Plus(Times(A, Xin, A), Times(minus1, A)),
                                    Omega
                                    )
                                )
                            )
                        )


class Example16():
    def __init__(self):

        # randomized matrix inversion 4

        # q << n
        n = 5000
        q = 500

        S = Matrix("S", (n, q))
        S.set_property(properties.FULL_RANK)
        S.set_property(properties.INPUT)

        A = Matrix("A", (n, n))
        A.set_property(properties.SYMMETRIC)
        A.set_property(properties.FULL_RANK)
        A.set_property(properties.INPUT)

        Xin = Matrix("Xin", (n, n))
        Xin.set_property(properties.FULL_RANK)
        Xin.set_property(properties.INPUT)

        Xout = Matrix("Xout", (n, n))
        Xout.set_property(properties.FULL_RANK)
        Xout.set_property(properties.OUTPUT)

        I = IdentityMatrix(n, n)
        minus1 = ConstantScalar(-1.0)

        subexpr = Times(S, Inverse(Times(Transpose(S), A, S)), Transpose(S))

        self.eqns = Equations(
                        Equal(Xout,
                            Plus(
                                subexpr,
                                Times(
                                    Plus(I, Times(minus1, subexpr, A)),
                                    Xin,
                                    Plus(I, Times(minus1, A, subexpr))
                                    )
                                )
                            )
                        )


class Example17():
    def __init__(self):

        # Stochastic Newton and Quasi-Newton for Large Linear Least-squares 1

        # l < n
        # n << m
        l = 625
        n = 1000
        m = 5000

        Wk = Matrix("Wk", (m, l))
        Wk.set_property(properties.FULL_RANK)
        Wk.set_property(properties.INPUT)

        A = Matrix("A", (m, n))
        A.set_property(properties.FULL_RANK)
        A.set_property(properties.INPUT)

        Bin = Matrix("Bin", (n, n))
        Bin.set_property(properties.SPD)
        Bin.set_property(properties.INPUT)

        Bout = Matrix("Bout", (n, n))
        Bout.set_property(properties.SPD)
        Bout.set_property(properties.OUTPUT)

        In = IdentityMatrix(n, n)
        Il = IdentityMatrix(l, l)
        k = Scalar("k")
        minus1 = ConstantScalar(-1.0)
        kminus1 = Plus(k, minus1)

        self.eqns = Equations(
                        Equal(Bout,
                            Times(
                                Times(k, Inverse(kminus1)),
                                Bin,
                                Plus(
                                    In,
                                    Times(
                                        minus1, Transpose(A), Wk,
                                        Inverse(Plus(
                                            Times(kminus1, Il),
                                            Times(Transpose(Wk), A, Bin, Transpose(A), Wk)
                                            )),
                                        Transpose(Wk), A, Bin
                                        )
                                    )
                                )
                            )
                        )


class Example18():
    def __init__(self):

        # Stochastic Newton and Quasi-Newton for Large Linear Least-squares 2

        # l < n
        # n << m
        l = 625
        n = 1000
        m = 5000

        Wk = Matrix("Wk", (m, l))
        Wk.set_property(properties.FULL_RANK)
        Wk.set_property(properties.INPUT)

        A = Matrix("A", (m, n))
        A.set_property(properties.FULL_RANK)
        A.set_property(properties.INPUT)

        Bout = Matrix("Bout", (n, n))
        Bout.set_property(properties.SPD)
        Bout.set_property(properties.OUTPUT)

        In = IdentityMatrix(n, n)
        Il = IdentityMatrix(l, l)
        lambda_ = Scalar("lambda")
        minus1 = ConstantScalar(-1.0)

        self.eqns = Equations(
                        Equal(Bout,
                            Times(
                                Times(ConstantScalar(1.0), Inverse(lambda_)),
                                Plus(
                                    In,
                                    Times(
                                        minus1, Transpose(A), Wk,
                                        Inverse(Plus(
                                            Times(lambda_, Il),
                                            Times(Transpose(Wk), A, Transpose(A), Wk)
                                            )),
                                        Transpose(Wk), A
                                        )
                                    )
                                )
                            )
                        )


class Example19():
    def __init__(self):

        # tikhonov regularization 1

        n = 1500
        m = 1000

        A = Matrix("A", (n, m), properties = [properties.INPUT, properties.FULL_RANK])
        I = IdentityMatrix(m, m)
        Gamma = Matrix("Gamma", (m , m), properties = [properties.INPUT, properties.FULL_RANK])
        b = Vector("b", (n, 1), properties = [properties.INPUT])
        x = Vector("x", (m, 1), properties = [properties.OUTPUT])

        self.init = lambda: derivation.special_properties.add_expression(Plus(Times(Transpose(A), A), Times(Transpose(Gamma), Gamma)), {properties.SPD})
        self.init()

        self.eqns = Equations(
                            Equal(x,
                                Times(Inverse(Plus(Times(Transpose(A), A), Times(Transpose(Gamma), Gamma))), Transpose(A), b)
                                )
                            )


class Example20():
    def __init__(self):

        # tikhonov regularization 2

        n = 1500
        m = 1000

        A = Matrix("A", (n, m), properties = [properties.INPUT, properties.FULL_RANK])
        I = IdentityMatrix(m, m)
        b = Vector("b", (n, 1), properties = [properties.INPUT])
        x = Vector("x", (m, 1), properties = [properties.OUTPUT])
        alpha = Scalar("alpha_sq")

        self.init = lambda: derivation.special_properties.add_expression(Plus(Times(Transpose(A), A), Times(alpha, I)), {properties.SPD})
        self.init()

        self.eqns = Equations(
                            Equal(x,
                                Times(Inverse(Plus(Times(Transpose(A), A), Times(alpha, I))), Transpose(A), b)
                                )
                            )


class Example21():
    def __init__(self):

        # tikhonov regularization 3

        n = 1500
        m = 1000

        A = Matrix("A", (n, m), properties = [properties.INPUT, properties.FULL_RANK])
        P = Matrix("P", (n, n), properties = [properties.INPUT, properties.FULL_RANK, properties.SYMMETRIC]) # covariance matrix
        Q = Matrix("Q", (m, m), properties = [properties.INPUT, properties.FULL_RANK, properties.SYMMETRIC]) # covariance matrix
        b = Vector("b", (n, 1), properties = [properties.INPUT])
        x0 = Vector("x0", (m, 1), properties = [properties.INPUT])
        x = Vector("x", (m, 1), properties = [properties.OUTPUT])

        self.eqns = Equations(
                            Equal(x,
                                Times(
                                    Inverse(Plus(Times(Transpose(A), P, A), Q)),
                                    Plus(Times(Transpose(A), P, b), Times(Q, x0))
                                    )
                                )
                            )


class Example22():
    def __init__(self):

        # tikhonov regularization 4

        n = 1500
        m = 1000

        A = Matrix("A", (n, m), properties = [properties.INPUT, properties.FULL_RANK])
        P = Matrix("P", (n, n), properties = [properties.INPUT, properties.FULL_RANK, properties.SYMMETRIC]) # covariance matrix
        Q = Matrix("Q", (m, m), properties = [properties.INPUT, properties.FULL_RANK, properties.SYMMETRIC]) # covariance matrix
        b = Vector("b", (n, 1), properties = [properties.INPUT])
        x0 = Vector("x0", (m, 1), properties = [properties.INPUT])
        x = Vector("x", (m, 1), properties = [properties.OUTPUT])
        minus1 = ConstantScalar(-1.0)

        self.eqns = Equations(
                            Equal(x,
                                Plus(x0,
                                    Times(
                                        Inverse(Plus(Times(Transpose(A), P, A), Q)),
                                        Transpose(A), P, Plus(b, Times(minus1, A, x0))
                                        )
                                    )
                                )
                            )


class Example23():
    def __init__(self):

        # linear MMSE estimator 1

        n = 1500
        m = 1000

        A = Matrix("A", (m, n), properties = [properties.INPUT, properties.FULL_RANK])
        Cx = Matrix("Cx", (n, n), properties = [properties.INPUT, properties.FULL_RANK, properties.SYMMETRIC]) # covariance matrix
        Cz = Matrix("Cz", (m, m), properties = [properties.INPUT, properties.FULL_RANK, properties.SYMMETRIC]) # covariance matrix
        y = Vector("y", (m, 1), properties = [properties.INPUT])
        x = Vector("x", (n, 1), properties = [properties.INPUT])
        xout = Vector("xout", (n, 1), properties = [properties.OUTPUT])
        minus1 = ConstantScalar(-1.0)

        self.eqns = Equations(
                        Equal(xout,
                            Plus(
                                Times(Cx, Transpose(A),
                                    Inverse(Plus(Times(A, Cx, Transpose(A)), Cz)),
                                    Plus(y, Times(minus1, A, x))
                                    ),
                                x)
                            )
                        )


class Example24():
    def __init__(self):

        # linear MMSE estimator 2

        n = 1500
        m = 1000

        A = Matrix("A", (m, n), properties = [properties.INPUT, properties.FULL_RANK])
        Cx = Matrix("Cx", (n, n), properties = [properties.INPUT, properties.FULL_RANK, properties.SYMMETRIC]) # covariance matrix
        Cz = Matrix("Cz", (m, m), properties = [properties.INPUT, properties.FULL_RANK, properties.SYMMETRIC]) # covariance matrix
        y = Vector("y", (m, 1), properties = [properties.INPUT])
        x = Vector("x", (n, 1), properties = [properties.INPUT])
        xout = Vector("xout", (n, 1), properties = [properties.OUTPUT])
        minus1 = ConstantScalar(-1.0)

        self.eqns = Equations(
                        Equal(xout,
                            Plus(
                                Times(
                                    Inverse(Plus(Times(Transpose(A), Inverse(Cz), A), Inverse(Cx))),
                                    Transpose(A), Inverse(Cz),
                                    Plus(y, Times(minus1, A, x))
                                    ),
                                x)
                            )
                        )


class Example25():
    def __init__(self):

        # linear MMSE estimator 3

        n = 1500
        m = 1000

        A = Matrix("A", (m, n), properties = [properties.INPUT, properties.FULL_RANK])
        Kin = Matrix("Kin", (n, m), properties = [properties.INPUT, properties.FULL_RANK])
        Cin = Matrix("Cin", (n, n), properties = [properties.INPUT, properties.FULL_RANK, properties.SYMMETRIC]) # covariance matrix
        Cz = Matrix("Cz", (m, m), properties = [properties.INPUT, properties.FULL_RANK, properties.SYMMETRIC]) # covariance matrix
        y = Vector("y", (m, 1), properties = [properties.INPUT])
        x = Vector("x", (n, 1), properties = [properties.INPUT])

        xout = Vector("xout", (n, 1), properties = [properties.OUTPUT])
        Cout = Matrix("Cout", (n, n), properties = [properties.OUTPUT])
        Kout = Matrix("Kout", (n, m), properties = [properties.OUTPUT])

        minus1 = ConstantScalar(-1.0)
        I = IdentityMatrix(n, n)

        self.eqns = Equations(
                        Equal(Kout, 
                            Times(Cin, Transpose(A), Inverse(Plus(Times(A, Cin, Transpose(A)), Cz)))),
                        Equal(xout,
                            Plus(x, Times(Kin, Plus(y, Times(minus1, A, x))))),
                        Equal(Cout, 
                            Times(Plus(I, Times(minus1, Kin, A)), Cin))
            )

class Example26():
    def __init__(self):

        # Kalman filter
        
        # TODO: R is symmetric positive semi-definite, which property?

        n = 1000
        m = 1000
        minus1 = ConstantScalar(-1.0)

        Kk_O = Matrix("Kk_O", (n, m), properties = [properties.FULL_RANK, properties.OUTPUT])
        Kk_I = Matrix("Kk_I", (n, m), properties = [properties.FULL_RANK, properties.INPUT])
        P_b = Matrix("P_b", (n, n), properties = [properties.INPUT, properties.SPD])
        P_a = Matrix("P_a", (n, n), properties = [properties.OUTPUT])
        H = Matrix("H", (m, n), properties = [properties.FULL_RANK, properties.INPUT])
        R = Matrix("R", (m, m), properties = [properties.FULL_RANK, properties.INPUT, properties.SYMMETRIC])
        I = IdentityMatrix(n, n)

        x_a = Vector("x_a", (n, 1), properties = [properties.OUTPUT])
        x_b = Vector("x_b", (n, 1), properties = [properties.INPUT])
        zk = Vector("zk", (m, 1), properties = [properties.INPUT])

        self.init = lambda: derivation.special_properties.add_expression(Plus(Times(H, P_b, Transpose(H)), R), {properties.SPD})
        self.init()

        self.eqns = Equations(
                            Equal(
                                Kk_O,
                                Times(
                                    P_b,
                                    Transpose(H),
                                    Inverse(Plus(Times(H, P_b, Transpose(H)), R))
                                )
                            ),
                            Equal(x_a, Plus(x_b, Times(Kk_I, Plus( zk, Times(minus1, H, x_b))))),
                            Equal(P_a, Times(Plus(I, Times(minus1, Kk_I, H)), P_b))
                        )