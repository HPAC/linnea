from ..algebra.expression import Symbol, Scalar, Vector, Matrix, ConstantScalar, \
                                 Equal, Plus, Times, Transpose, Inverse, \
                                 InverseTranspose, InverseConjugate, \
                                 InverseConjugateTranspose, \
                                 ConjugateTranspose, Index, IdentityMatrix

from ..algebra.properties import Property

from ..algebra.equations import Equations

from .. import derivation

from ..derivation import special_properties


class Example01():
    def __init__(self):

        # least squares

        n = 2500
        m = 500

        X = Matrix("X", (n, m))
        X.set_property(Property.FULL_RANK)

        y = Vector("y", (n, 1))

        b = Vector("b", (m, 1))

        # b = (X^T X)^-1 X^T y
        self.eqns = Equations(Equal(b, Times(Inverse(Times(Transpose(X), X)), Transpose(X), y)))


class Example02():
    def __init__(self):

        # generalized least squares

        n = 2500
        m = 500

        S = Matrix("S", (n, n))
        S.set_property(Property.SPD)

        X = Matrix("X", (n, m))
        X.set_property(Property.FULL_RANK)

        z = Vector("z", (m, 1))

        y = Vector("y", (n, 1))

        self.eqns = Equations(Equal(z, Times(Inverse(Times(Transpose(X), Inverse(S), X ) ), Transpose(X), Inverse(S), y)))


class Example03():
    def __init__(self, ):

        # optimization problem 1

        n = 2000
        m = 1000

        A = Matrix("A", (m, n))
        A.set_property(Property.FULL_RANK)

        W = Matrix("W", (n, n))
        # W is positive
        W.set_property(Property.FULL_RANK)
        W.set_property(Property.DIAGONAL)
        W.set_property(Property.SPD)

        b = Vector("b", (m, 1))

        c = Vector("c", (n, 1))

        x = Vector("x", (n, 1))

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

        n = 2000
        m = 1000

        A = Matrix("A", (m, n))
        A.set_property(Property.FULL_RANK)

        W = Matrix("W", (n, n))
        # W is positive
        W.set_property(Property.FULL_RANK)
        W.set_property(Property.DIAGONAL)
        W.set_property(Property.SPD)

        b = Vector("b", (m, 1))

        c = Vector("c", (n, 1))

        x = Vector("x", (n, 1))

        xf = Vector("xf", (n, 1))

        xo = Vector("xo", (n, 1))

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
    def __init__(self):

        N = 2000

        # signal processing

        A = Matrix("A", (N, N), properties = [Property.FULL_RANK])
        # A - tridiagonal, full rank, something even more specific
        B = Matrix("B", (N, N), properties = [Property.FULL_RANK])
        # B - tridiagonal, full rank, something even more specific
        R = Matrix("R", (N - 1, N), properties = [Property.FULL_RANK, Property.UPPER_TRIANGULAR])
        # R - upper bidiagonal
        L = Matrix("L", (N - 1, N - 1), properties = [Property.FULL_RANK, Property.DIAGONAL])

        y = Vector("y", (N, 1))
        x = Vector("x", (N, 1))

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

        n = 2000
        m = 200
        k = 2000

        L00 = Matrix("L00", (n, n), properties = [Property.FULL_RANK, Property.LOWER_TRIANGULAR])
        L11 = Matrix("L11", (m, m), properties = [Property.FULL_RANK, Property.LOWER_TRIANGULAR])
        L22 = Matrix("L22", (k, k), properties = [Property.FULL_RANK, Property.LOWER_TRIANGULAR])
        L21 = Matrix("L21", (k, m), properties = [Property.FULL_RANK])
        L10 = Matrix("L10", (m, n), properties = [Property.FULL_RANK])
        L20 = Matrix("L20", (k, n), properties = [Property.FULL_RANK])

        X21 = Matrix("X21", (k, m))
        X11 = Matrix("X11", (m, m))
        X10 = Matrix("X10", (m, n))
        
        X20 = Matrix("X20", (k, n))
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

        N = 200
        msd = 2000 # p*nsd, p < 1 (percentage)
        nsd = 2000

        # TODO there is a dimension mismatch here for msd != nsd

        minus1 = ConstantScalar(-1.0)

        B = Matrix("B", (nsd, nsd), properties = [Property.SPSD]) # covariance matrix
        H = Matrix("H", (msd, nsd), properties = [Property.FULL_RANK])
        R = Matrix("R", (msd, msd), properties = [Property.SPSD]) # covariance matrix
        Y = Matrix("Y", (msd, N), properties = [Property.FULL_RANK])
        Xb = Matrix("Xb", (nsd, N), properties = [Property.FULL_RANK])
        Xa = Matrix("Xa", (nsd, N), properties = [Property.FULL_RANK])

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

        N = 200
        m = 1000 # p*nsd, p < 1 (percentage)
        n = 5000

        minus1 = ConstantScalar(-1.0)

        B = Matrix("B", (n, n), properties = [Property.SPSD]) # covariance matrix
        H = Matrix("H", (m, n), properties = [Property.FULL_RANK])
        R = Matrix("R", (m, m), properties = [Property.SPSD]) # covariance matrix
        Y = Matrix("Y", (m, N), properties = [Property.FULL_RANK])
        Xb = Matrix("Xb", (n, N), properties = [Property.FULL_RANK])
        dX = Matrix("dX", (n, N), properties = [Property.FULL_RANK])

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

        N = 200
        m = 1000 # p*nsd, p < 1 (percentage)
        n = 5000

        minus1 = ConstantScalar(-1.0)

        H = Matrix("H", (m, n), properties = [Property.FULL_RANK])
        R = Matrix("R", (m, m), properties = [Property.SPSD]) # covariance matrix
        Y = Matrix("Y", (m, N), properties = [Property.FULL_RANK])
        X = Matrix("X", (n, n), properties = [Property.FULL_RANK, Property.LOWER_TRIANGULAR])
        Xb = Matrix("Xb", (n, N), properties = [Property.FULL_RANK])
        dX = Matrix("dX", (n, N), properties = [Property.FULL_RANK])

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

        n = 5000
        m = 1000

        minus1 = ConstantScalar(-1.0)
        lambda_ = Scalar("lambda")
        lambda_.set_property(Property.POSITIVE)
        sigma_ = Scalar("sigma_sq")
        sigma_.set_property(Property.POSITIVE)

        H = Matrix("H", (m, n), properties = [Property.FULL_RANK])
        I = IdentityMatrix(n, n)

        v_k = Vector("v_k", (n, 1))
        u_k = Vector("u_k", (n, 1))
        y = Vector("y", (m, 1))
        x = Vector("x", (n, 1))


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

        n = 5000
        m = 1000

        minus1 = ConstantScalar(-1.0)
        lambda_ = Scalar("lambda")
        lambda_.set_property(Property.POSITIVE)
        sigma_ = Scalar("sigma_sq")
        sigma_.set_property(Property.POSITIVE)

        H = Matrix("H", (m, n), properties = [Property.FULL_RANK])
        H_dag = Matrix("H_dag", (n, m), properties = [Property.FULL_RANK])
        I = IdentityMatrix(n, n)

        y_k = Vector("y_k", (n, 1))
        y = Vector("y", (m, 1))
        x = Vector("x", (n, 1))

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
                            H_dag,
                            h_dag
                        ),
                        Equal(
                            y_k,
                            Plus(
                                Times(
                                    H_dag,
                                    y
                                ),
                                Times(
                                    Plus(
                                        I,
                                        Times(
                                            minus1,
                                            H_dag,
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

        # randomized matrix inversion 1

        # q << n
        n = 5000
        q = 500

        W = Matrix("W", (n, n))
        W.set_property(Property.SPD)

        S = Matrix("S", (n, q))
        S.set_property(Property.FULL_RANK)

        A = Matrix("A", (n, n))
        A.set_property(Property.FULL_RANK)

        Lambda = Matrix("Lambda", (n, n))
        Lambda.set_property(Property.SYMMETRIC)
        Lambda.set_property(Property.FULL_RANK)

        Xin = Matrix("Xin", (n, n))
        Xin.set_property(Property.FULL_RANK)

        Xout = Matrix("Xout", (n, n))

        I = IdentityMatrix(n, n)
        minus1 = ConstantScalar(-1.0)

        self.eqns = Equations(
                        Equal(Lambda,
                                Times(S, Inverse(Times(Transpose(S), A, W, Transpose(A), S)), Transpose(S))),
                        Equal(Xout,
                            Plus(Xin,
                                Times(
                                    W, Transpose(A), Lambda,
                                    Plus(I, Times(minus1, A, Xin))
                                    )
                                )
                            )
                        )

class Example13():
    def __init__(self):

        # randomized matrix inversion 2

        # q << n
        n = 5000
        q = 500

        W = Matrix("W", (n, n))
        W.set_property(Property.SPD)

        S = Matrix("S", (n, q))
        S.set_property(Property.FULL_RANK)

        A = Matrix("A", (n, n))
        A.set_property(Property.FULL_RANK)

        Lambda = Matrix("Lambda", (n, n))
        Lambda.set_property(Property.SYMMETRIC)
        Lambda.set_property(Property.FULL_RANK)

        Xin = Matrix("Xin", (n, n))
        Xin.set_property(Property.FULL_RANK)

        Xout = Matrix("Xout", (n, n))

        I = IdentityMatrix(n, n)
        minus1 = ConstantScalar(-1.0)

        self.eqns = Equations(
                        Equal(Lambda,
                            Times(S, Inverse(Times(Transpose(S), Transpose(A), W, A, S)), Transpose(S))),
                        Equal(Xout,
                            Plus(Xin,
                                Times(
                                    Plus(I, Times(minus1, Xin, Transpose(A))),
                                    Lambda, Transpose(A), W
                                    )
                                )
                            )
                        )


class Example14():
    def __init__(self):

        # randomized matrix inversion 3

        # q << n
        n = 5000
        q = 500

        W = Matrix("W", (n, n))
        W.set_property(Property.SPD)

        S = Matrix("S", (n, q))
        S.set_property(Property.FULL_RANK)

        A = Matrix("A", (n, n))
        A.set_property(Property.SYMMETRIC)
        A.set_property(Property.FULL_RANK)

        Lambda = Matrix("Lambda", (n, n))
        Lambda.set_property(Property.SYMMETRIC)
        Lambda.set_property(Property.FULL_RANK)

        Theta = Matrix("Theta", (n, n))
        Theta.set_property(Property.FULL_RANK)

        Mk = Matrix("Mk", (n, n))
        Mk.set_property(Property.FULL_RANK)

        Xin = Matrix("Xin", (n, n))
        Xin.set_property(Property.SYMMETRIC)
        Xin.set_property(Property.FULL_RANK)

        Xout = Matrix("Xout", (n, n))

        I = IdentityMatrix(n, n)
        minus1 = ConstantScalar(-1.0)


        self.eqns = Equations(
                        Equal(Lambda, Times(S, Inverse(Times(Transpose(S), A, W, A, S)), Transpose(S))),
                        Equal(Theta, Times(Lambda, A, W)),
                        Equal(Mk, Plus(Times(Xin, A), Times(minus1, I))),
                        Equal(Xout,
                            Plus(Xin,
                                Times(minus1, Mk, Theta),
                                Times(minus1, Transpose(Times(Mk, Theta))),
                                Times(
                                    Transpose(Theta),
                                    Plus(Times(A, Xin, A), Times(minus1, A)),
                                    Theta
                                    )
                                )
                            )
                        )


class Example15():
    def __init__(self):

        # randomized matrix inversion 4

        # q << n
        n = 5000
        q = 500

        S = Matrix("S", (n, q))
        S.set_property(Property.FULL_RANK)

        A = Matrix("A", (n, n))
        A.set_property(Property.SPD)
        A.set_property(Property.FULL_RANK)

        Xin = Matrix("Xin", (n, n))
        Xin.set_property(Property.SYMMETRIC)
        Xin.set_property(Property.FULL_RANK)

        Xout = Matrix("Xout", (n, n))
        Xout.set_property(Property.FULL_RANK)

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


class Example16():
    def __init__(self):

        # Stochastic Newton and Quasi-Newton for Large Linear Least-squares 1

        # l < n
        # n << m
        l = 625
        n = 1000
        m = 5000

        Wk = Matrix("Wk", (m, l))
        Wk.set_property(Property.FULL_RANK)

        A = Matrix("A", (m, n))
        A.set_property(Property.FULL_RANK)

        Bin = Matrix("Bin", (n, n))
        Bin.set_property(Property.SPD)

        Bout = Matrix("Bout", (n, n))
        Bout.set_property(Property.SPD)

        In = IdentityMatrix(n, n)
        Il = IdentityMatrix(l, l)
        k = Scalar("k")
        k.set_property(Property.POSITIVE)
        minus1 = ConstantScalar(-1.0)
        kminus1 = Plus(k, minus1) # this is positive too, but using special_properties doesn't work

        # This works:
        # special_properties.add_expression(Plus(Times(kminus1, Il), Times(Transpose(Wk), A, Bin, Transpose(A), Wk)), [Property.SPD])
        
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


class Example17():
    def __init__(self):

        # Stochastic Newton and Quasi-Newton for Large Linear Least-squares 2

        # l < n
        # n << m
        l = 625
        n = 1000
        m = 5000

        Wk = Matrix("Wk", (m, l))
        Wk.set_property(Property.FULL_RANK)

        A = Matrix("A", (m, n))
        A.set_property(Property.FULL_RANK)

        B = Matrix("B", (n, n))
        B.set_property(Property.SPD)

        In = IdentityMatrix(n, n)
        Il = IdentityMatrix(l, l)
        lambda_ = Scalar("lambda")
        lambda_.set_property(Property.POSITIVE)
        minus1 = ConstantScalar(-1.0)

        self.eqns = Equations(
                        Equal(B,
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


class Example18():
    def __init__(self):

        # tikhonov regularization 1

        n = 3000
        m = 200

        A = Matrix("A", (n, m), properties = [Property.FULL_RANK])
        Gamma = Matrix("Gamma", (m , m), properties = [Property.FULL_RANK])
        b = Vector("b", (n, 1))
        x = Vector("x", (m, 1))

        self.eqns = Equations(
                            Equal(x,
                                Times(Inverse(Plus(Times(Transpose(A), A), Times(Transpose(Gamma), Gamma))), Transpose(A), b)
                                )
                            )


class Example19():
    def __init__(self):

        # tikhonov regularization 2

        n = 3000
        m = 200

        A = Matrix("A", (n, m), properties = [Property.FULL_RANK])
        I = IdentityMatrix(m, m)
        b = Vector("b", (n, 1))
        x = Vector("x", (m, 1))
        alpha = Scalar("alpha_sq")
        alpha.set_property(Property.POSITIVE)

        self.eqns = Equations(
                            Equal(x,
                                Times(Inverse(Plus(Times(Transpose(A), A), Times(alpha, I))), Transpose(A), b)
                                )
                            )


class Example20():
    def __init__(self):

        # tikhonov regularization 3

        n = 3000
        m = 200

        A = Matrix("A", (n, m), properties = [Property.FULL_RANK])
        P = Matrix("P", (n, n), properties = [Property.SPD]) # covariance matrix
        Q = Matrix("Q", (m, m), properties = [Property.SPD]) # covariance matrix
        b = Vector("b", (n, 1))
        x0 = Vector("x0", (m, 1))
        x = Vector("x", (m, 1))

        self.eqns = Equations(
                            Equal(x,
                                Times(
                                    Inverse(Plus(Times(Transpose(A), P, A), Q)),
                                    Plus(Times(Transpose(A), P, b), Times(Q, x0))
                                    )
                                )
                            )


class Example21():
    def __init__(self):

        # tikhonov regularization 4

        n = 3000
        m = 200

        A = Matrix("A", (n, m), properties = [Property.FULL_RANK])
        P = Matrix("P", (n, n), properties = [Property.SPD]) # covariance matrix
        Q = Matrix("Q", (m, m), properties = [Property.SPD]) # covariance matrix
        b = Vector("b", (n, 1))
        x0 = Vector("x0", (m, 1))
        x = Vector("x", (m, 1))
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


class Example22():
    def __init__(self):

        # linear MMSE estimator 1

        n = 2000
        m = 1500

        A = Matrix("A", (m, n), properties = [Property.FULL_RANK])
        Cx = Matrix("Cx", (n, n), properties = [Property.SPD]) # covariance matrix
        Cz = Matrix("Cz", (m, m), properties = [Property.SPD]) # covariance matrix
        y = Vector("y", (m, 1))
        x = Vector("x", (n, 1))
        xout = Vector("xout", (n, 1))
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


class Example23():
    def __init__(self):

        # linear MMSE estimator 2

        n = 2000
        m = 1500

        A = Matrix("A", (m, n), properties = [Property.FULL_RANK])
        Cx = Matrix("Cx", (n, n), properties = [Property.SPSD]) # covariance matrix
        Cz = Matrix("Cz", (m, m), properties = [Property.SPSD]) # covariance matrix
        y = Vector("y", (m, 1))
        x = Vector("x", (n, 1))
        xout = Vector("xout", (n, 1))
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


class Example24():
    def __init__(self):

        # linear MMSE estimator 3

        n = 2000
        m = 1500

        A = Matrix("A", (m, n), properties = [Property.FULL_RANK])
        K = Matrix("K", (n, m), properties = [Property.FULL_RANK])
        Cin = Matrix("Cin", (n, n), properties = [Property.SPD]) # covariance matrix
        Cz = Matrix("Cz", (m, m), properties = [Property.SPD]) # covariance matrix
        y = Vector("y", (m, 1))
        x = Vector("x", (n, 1))

        xout = Vector("xout", (n, 1))
        Cout = Matrix("Cout", (n, n)) # covariance matrix

        minus1 = ConstantScalar(-1.0)
        I = IdentityMatrix(n, n)

        self.eqns = Equations(
                        Equal(K, 
                            Times(Cin, Transpose(A), Inverse(Plus(Times(A, Cin, Transpose(A)), Cz)))),
                        Equal(xout,
                            Plus(x, Times(K, Plus(y, Times(minus1, A, x))))),
                        Equal(Cout, 
                            Times(Plus(I, Times(minus1, K, A)), Cin))
            )

class Example25():
    def __init__(self):

        # Kalman filter

        n = 400
        m = 400
        minus1 = ConstantScalar(-1.0)

        Kk = Matrix("Kk", (n, m), properties = [Property.FULL_RANK])
        P_b = Matrix("P_b", (n, n), properties = [Property.SPD])
        P_a = Matrix("P_a", (n, n))
        H = Matrix("H", (m, n), properties = [Property.FULL_RANK])
        R = Matrix("R", (m, m), properties = [Property.SPSD]) # covariance matrix
        I = IdentityMatrix(n, n)

        x_a = Vector("x_a", (n, 1))
        x_b = Vector("x_b", (n, 1))
        zk = Vector("zk", (m, 1))

        self.eqns = Equations(
                            Equal(
                                Kk,
                                Times(
                                    P_b,
                                    Transpose(H),
                                    Inverse(Plus(Times(H, P_b, Transpose(H)), R))
                                )
                            ),
                            Equal(x_a, Plus(x_b, Times(Kk, Plus( zk, Times(minus1, H, x_b))))),
                            Equal(P_a, Times(Plus(I, Times(minus1, Kk, H)), P_b))
                        )