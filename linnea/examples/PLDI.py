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

        # 1.1
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

        # 1.2
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

        # 1.3
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

        # 1.4
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

        # 1.5
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

        # 1.6
        # inversion of lower triangular matrix

        n = 1000
        m = 1000
        k = 1000

        L00 = Matrix("L00", (n, n), properties = [properties.FULL_RANK, properties.LOWER_TRIANGULAR, properties.INPUT])
        L11 = Matrix("L11", (m, m), properties = [properties.FULL_RANK, properties.LOWER_TRIANGULAR, properties.INPUT])
        L22 = Matrix("L22", (k, k), properties = [properties.FULL_RANK, properties.LOWER_TRIANGULAR, properties.INPUT])
        L21 = Matrix("L21", (k, m), properties = [properties.FULL_RANK, properties.INPUT])
        L10 = Matrix("L10", (m, n), properties = [properties.FULL_RANK, properties.INPUT])

        X21 = Matrix("X21", (k, m), properties = [properties.OUTPUT])
        X11 = Matrix("X11", (m, m), properties = [properties.OUTPUT])
        X10Input = Matrix("X10Input", (m, n), properties = [properties.INPUT])
        X10Output = Matrix("X10Output", (m, n), properties = [properties.OUTPUT])
        X20Input = Matrix("X20Input", (k, n), properties = [properties.INPUT])
        X20Output = Matrix("X20Output", (k, n), properties = [properties.OUTPUT])
        minus1 = ConstantScalar(-1.0)

        self.eqns = Equations(
                            Equal(
                                X10Output,
                                Times(
                                    X10Input,
                                    Inverse(
                                        L00
                                    )
                                )
                            ),
                            Equal(
                                X20Output,
                                Plus(
                                    X20Input,
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

        # 1.7
        # local assimilation for parallel ensemble Kalman filter based on modified Cholesky decomposition

        N = 1000
        msd = 1000
        nsd = 1000

        minus1 = ConstantScalar(-1)

        B = Matrix("B", (N, N), properties = [properties.FULL_RANK, properties.SYMMETRIC, properties.INPUT])
        # B is a covariance matrix (symmetric positive semi-definite)
        H = Matrix("H", (msd, N), properties = [properties.FULL_RANK, properties.INPUT])
        R = Matrix("R", (msd, msd), properties = [properties.FULL_RANK, properties.SYMMETRIC, properties.INPUT])
        # R is a covariance matrix (symmetric positive semi-definite)
        Y = Matrix("Y", (msd, N), properties = [properties.FULL_RANK, properties.INPUT])
        Xb = Matrix("Xb", (nsd, N), properties = [properties.FULL_RANK, properties.INPUT])
        Xa = Matrix("X", (nsd, N), properties = [properties.FULL_RANK, properties.OUTPUT])

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