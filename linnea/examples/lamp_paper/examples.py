from ...algebra.expression import Symbol, Scalar, Vector, Matrix, ConstantScalar, \
                                Equal, Plus, Times, Transpose, Inverse, \
                                InverseTranspose, InverseConjugate, \
                                InverseConjugateTranspose, \
                                ConjugateTranspose, Index, IdentityMatrix

from ...algebra.properties import Property as properties

from ...algebra.equations import Equations

from ... import derivation

sONE = 1
sZERO = 0

class LeastSquares_7_1_1(object):
    def __init__(self, n = 1500, m = 1000):

        # TAGS
        # least squares
        # problem 7.1.1 in the paper
        # equation: (X^T X)^-1 X^T y
        # where X is of size (n, m) and has full rank
        # y is of size (n, 1)
        # and n > m
        # vector solution b is of size (m, 1)
        # Example 007 in examples.py


        X = Matrix("X", (n, m), properties = [properties.INPUT, properties.FULL_RANK])

        y = Vector("y", (n, sONE), properties = [properties.INPUT])

        b = Vector("b", (m, sONE), properties = [properties.OUTPUT])

        # b = (X^T X)^-1 X^T y
        self.eqns = Equations(Equal(b, Times(Inverse(Times(Transpose(X), X)), Transpose(X), y)))

class LMMSE_7_1_2(object):
    def __init__(self, l = 1000, m = 500, n = 1000):

        # TAGS
        # mean square error
        # chen2017


        H = Matrix("H", (l + m, n), properties = [properties.INPUT])
        I = IdentityMatrix(n, n)
        y = Vector("y", (l + m, sONE), properties = [properties.INPUT])
        x = Vector("x", (n, sONE), properties = [properties.OUTPUT])
        sigma = Scalar("sigma")

        # x = (H^H*H + \sigma^2 * I)^-1 * H^H * y
        self.eqns = Equations(
                        Equal(x,
                            Times(
                                Inverse(
                                    Plus(
                                        Times(
                                            ConjugateTranspose(H),
                                            H
                                        ),
                                        Times(
                                            sigma,
                                            I
                                        )
                                    )
                                ),
                                ConjugateTranspose(H),
                                y
                            )
                        )
                    )

class Generalized_LeastSquares_7_1_3(object):
    def __init__(self, n = 1500, m = 1000):

        # TAGS
        # least squares
        # problem 7.1.3 in the paper
        # equation: b = (X^T M^-1 X)^-1 X^T M^-1 y
        # where X is of size (n, m) and has full rank
        # M is of size (m,m) and is SPD
        # y is of size (n, 1)
        # and n > m
        # vector solution b is of size (m, 1)
        # Example 007 in examples.py

        M = Matrix("M", (n, n), properties = [properties.SPD, properties.INPUT])

        X = Matrix("X", (n, m), properties = [properties.FULL_RANK, properties.INPUT])
        y = Matrix("y", (n, sONE), properties = [properties.INPUT])

        b = Vector("b", (m, sONE), properties = [properties.OUTPUT])

        self.eqns = Equations(Equal(b, Times(Inverse(Times(Transpose(X), Inverse(M), X ) ), Transpose(X), Inverse(M), y)))

class Optimization_Problem_7_1_4(object):
    def __init__(self, n = 1500, m = 1000):

        # TAGS
        # optimization problem
        # straszak2015
        # problem 7.1.4 in the paper
        # equation: x = W(A^T(AWA^T)^-1*b - c)
        # where A is of size (m, n) and has full rank
        # W is of size (n, n) and is diagonal, full rank, positive
        # b is of size (m, 1)
        # c is of size (n, 1)
        # and vector solution x is of size (n, 1)
        # and n > m
        # Example 126 in examples.py

        # n > m

        A = Matrix("A", (m, n), properties = [properties.FULL_RANK, properties.INPUT])

        W = Matrix("W", (n, n), properties = [properties.DIAGONAL, properties.FULL_RANK, properties.SPD, properties.INPUT])

        # W is positive
        #W.set_property(properties.FULL_RANK)
        #W.set_property(properties.DIAGONAL)
        #W.set_property(properties.SPD)
        #W.set_property(properties.SYMMETRIC)
        #W.set_property(properties.NON_SINGULAR)
        #W.set_property(properties.INPUT)

        b = Vector("b", (m, sONE), properties = [properties.INPUT])
        c = Vector("c", (n, sONE), properties = [properties.INPUT])
        x = Vector("x", (n, sONE), properties = [properties.OUTPUT])

        minusone = ConstantScalar(-1.0)


        self.eqns = Equations(
                            Equal(
                                x,
                                Times(
                                    W,
                                    Plus(
                                        Times(Transpose(A), Inverse(Times(A, W, Transpose(A))), b),
                                        Times(minusone, c)
                                        )
                                    )
                                )
                            )

class Signal_Processing_7_1_5(object):
    def __init__(self, N = 1000):
        # TAGS
        # signal processing
        # ding2016
        # problem 7.1.5 in the paper
        # equation: x = (A^-T B^T BA^-1 + R^TDR)^-1 A^-T B^T BA^-1 y
        # A, B - tridiagonal, full rank
        # R - upper bidiagonal
        # D - diagonal

        A = Matrix("A", (N, N), properties = [properties.FULL_RANK, properties.INPUT])
        B = Matrix("B", (N, N), properties = [properties.FULL_RANK, properties.INPUT])
        R = Matrix("R", (N - 1, N), properties = [properties.UPPER_TRIANGULAR, properties.INPUT])
        L = Matrix("L", (N - 1, N - 1), properties = [properties.DIAGONAL, properties.INPUT])

        y = Vector("y", (N, sONE), properties = [properties.INPUT])
        x = Vector("x", (N, sONE), properties = [properties.OUTPUT])

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

class Lower_Triangular_Inversion_7_1_6(object):
    def __init__(self, n = 1000, m = 1000, k = 1000):
        # TAGS
        # signal processing
        # problem 7.1.6 in the paper
        # equation: x = (A^-T B^T BA^-1 + R^TDR)^-1 A^-T B^T BA^-1 y
        # input-output variables are modeled as two seperate variables

        L00  = Matrix("L00", (n, n), properties = [properties.LOWER_TRIANGULAR, properties.INPUT])
        L11  = Matrix("L11", (m, m), properties = [properties.LOWER_TRIANGULAR, properties.INPUT])
        L22  = Matrix("L22", (k, k), properties = [properties.LOWER_TRIANGULAR, properties.INPUT])
        L21  = Matrix("L21", (k, m), properties = [properties.INPUT])
        L10  = Matrix("L10", (m, n), properties = [properties.INPUT])

        X21  = Matrix("X21", (k, m), properties = [properties.OUTPUT])
        X11  = Matrix("X11", (m, m), properties = [properties.OUTPUT])
        X10Input  = Matrix("X10Input", (m, n), properties = [properties.INPUT])
        X10Output  = Matrix("X10Output", (m, n), properties = [properties.OUTPUT])
        X20Input  = Matrix("X20Input", (k, n), properties = [properties.INPUT])
        X20Output  = Matrix("X20Output", (k, n), properties = [properties.OUTPUT])
        minusone = ConstantScalar(-1.0)

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
                                    minusone,
                                    Inverse(L22),
                                    L21
                                )
                            )
                    )

#FIXME: bug
class Local_Assimilation_Kalmar_7_1_7(object):
    def __init__(self, N = 1000, msd = 1000, nsd = 1000):
        # nino2016
        # something something Kalmar filter
        # problem 7.1.7 in the paper
        # equation: Xa = Xb+ (B^-1 + H^T*R^-1*H)^-1(Y-H*Xb)

        minusone = ConstantScalar(-1)

        B = Matrix("B", (N, N), properties = [properties.INPUT])
        H = Matrix("H", (msd, N), properties = [properties.INPUT])
        R = Matrix("R", (msd, msd), properties = [properties.INPUT])
        Y = Matrix("Y", (msd, N), properties = [properties.INPUT])
        Xb = Matrix("Xb", (nsd, N), properties = [properties.INPUT])
        Xa = Matrix("X", (nsd, N), properties = [properties.OUTPUT])

        self.eqns = Equations(
                            Equal(
                                Xa,
                                Plus(
                                    Xb,
                                    Times(
                                        Inverse(Plus(Inverse(B), Times(Transpose(H), Inverse(R), H))),
                                        # Times(H, Xb, minusone)
                                        Plus(Y, Times(H, Xb))
                                        )
                                    )
                                )
                            )

class Common_Subexpr_7_2_1(object):
    def __init__(self, n=1000):
        # similar to example 124
        # equation: X = ABCD, Y = BCDE, Z = BCDEF


        A = Matrix("A", (n, n), properties = [properties.INPUT])
        B = Matrix("B", (n, n), properties = [properties.INPUT])
        C = Matrix("C", (n, n), properties = [properties.INPUT])
        D = Matrix("D", (n, n), properties = [properties.INPUT])
        E = Matrix("E", (n, n), properties = [properties.INPUT])
        F = Matrix("F", (n, n), properties = [properties.INPUT])

        X = Matrix("X", (n, n), properties = [properties.OUTPUT])
        Y = Matrix("Y", (n, n), properties = [properties.OUTPUT])
        Z = Matrix("Z", (n, n), properties = [properties.OUTPUT])

        self.eqns = Equations(
                            Equal(X, Times(A, B, C, D)),
                            Equal(Y, Times(B, C, D, E)),
                            Equal(Z, Times(B, C, D, E, F)),
                            )

class Common_Subexpr_7_2_2(object):
    def __init__(self, n = 1500, m = 100):
        # equation: X = AB*AB
        # optimal solution A((BA)B when n >> m


        A = Matrix("A", (n, m), properties = [properties.INPUT])
        B = Matrix("B", (m, n), properties = [properties.INPUT])

        X = Matrix("X", (n, n), properties = [properties.OUTPUT])

        self.eqns = Equations(
                            Equal(X, Times(A, B, A, B)),
                            )

class Common_Subexpr_7_2_3(object):
    def __init__(self, n = 1500):
        # equation: X = AS, Y = SA^T
        # optimal solution uses AS = (SA^T)^T


        A = Matrix("A", (n, n), properties = [properties.INPUT])
        S = Matrix("S", (n, n), properties = [properties.INPUT, properties.SYMMETRIC])

        X = Matrix("X", (n, n), properties = [properties.OUTPUT])
        Y = Matrix("Y", (n, n), properties = [properties.OUTPUT])

        self.eqns = Equations(
                            Equal(X, Times(A, S)),
                            Equal(Y, Times(S, Transpose(A)))
                            )

class Overlap_Common_Subexpr_7_2_4(object):
    def __init__(self, m = 1500, n = 100):
        # equation: X = AB, Y = ABC, Z = BC
        # m >> n
        # optimal solution reuses BC


        A = Matrix("A", (n, n), properties = [properties.INPUT])
        B = Matrix("B", (n, n), properties = [properties.INPUT])
        C = Matrix("C", (n, m), properties = [properties.INPUT])

        X = Matrix("X", (n, n), properties = [properties.OUTPUT])
        Y = Matrix("Y", (n, m), properties = [properties.OUTPUT])
        Z = Matrix("Z", (n, m), properties = [properties.OUTPUT])

        self.eqns = Equations(
                            Equal(X, Times(A, B)),
                            Equal(Y, Times(A, B, C)),
                            Equal(Z, Times(B, C))
                            )

class Rewrite_Distributivity_Base(object):
    def create_matrices(self, n):
        A = Matrix("A", (n, n), properties = [properties.INPUT])
        B = Matrix("B", (n, n), properties = [properties.INPUT])
        C = Matrix("C", (n, n), properties = [properties.INPUT])
        D = Matrix("D", (n, n), properties = [properties.INPUT])
        E = Matrix("E", (n, n), properties = [properties.INPUT])
        F = Matrix("F", (n, n), properties = [properties.INPUT])

        X = Matrix("X", (n, n), properties = [properties.OUTPUT])
        return [A, B, C, D, E, F, X]

class Rewrite_Distributivity_7_2_5_1(Rewrite_Distributivity_Base):
    def __init__(self, n = 1500):
        # compute independently to test assignment 
        # X1 = AE + (A + B)(C + D)
        [A, B, C, D, E, F, X] = self.create_matrices(n)
        self.eqns = Equations(
                            Equal(X, Plus(
                                Times(
                                    A,
                                    E
                                ),
                                Times(
                                    Plus(A, B),
                                    Plus(C, D)
                                )
                            ))
                    )

class Rewrite_Distributivity_7_2_5_2(Rewrite_Distributivity_Base):
    def __init__(self, n = 1500):
        # X2 = AE + A(C + D) + B(C + D)
        [A, B, C, D, E, F, X] = self.create_matrices(n)
        self.eqns = Equations(
                            Equal(X, Plus(
                                Times(
                                    A,
                                    E
                                ),
                                Times(
                                    A,
                                    Plus(C, D)
                                ),
                                Times(
                                    B,
                                    Plus(C, D)
                                )
                            ))
                    )

class Rewrite_Distributivity_7_2_5_3(Rewrite_Distributivity_Base):
    def __init__(self, n = 1500):
        # X3 = A(C + D + E) + B(C + D)
        [A, B, C, D, E, F, X] = self.create_matrices(n)
        self.eqns = Equations(
                            Equal(X, Plus(
                                Times(
                                    A,
                                    Plus(C, D, E),
                                ),
                                Times(
                                    B,
                                    Plus(C, D)
                                )
                            ))
                    )

class Rewrite_Distributivity_7_2_5_4(Rewrite_Distributivity_Base):
    def __init__(self, n = 1500):
        # X4 = AC + AD + AE + BC + BD
        [A, B, C, D, E, F, X] = self.create_matrices(n)
        self.eqns = Equations(
                            Equal(X, Plus(
                                Times(A, C),
                                Times(A, D),
                                Times(A, E),
                                Times(B, C),
                                Times(B, D),
                            ))
                    )

class Matrix_Chain_7_2_6(object):
    def __init__(self, n1 = 600, n2 = 200, n3 = 1500, n4 = 800, n5 = 1000):
        # equation: W = ABCDEF
        # A is lower triangular, 
        # D, E are upper triangular

        A = Matrix("A", (n1, n1), properties = [properties.INPUT, properties.LOWER_TRIANGULAR])
        B = Matrix("B", (n1, n2), properties = [properties.INPUT])
        C = Matrix("C", (n2, n3), properties = [properties.INPUT])
        D = Matrix("D", (n3, n3), properties = [properties.INPUT, properties.UPPER_TRIANGULAR])
        E = Matrix("E", (n3, n4), properties = [properties.INPUT, properties.UPPER_TRIANGULAR])
        F = Matrix("F", (n4, n5), properties = [properties.INPUT])

        W = Matrix("W", (n1, n5), properties = [properties.OUTPUT])

        self.eqns = Equations(
                            Equal(W, Times(A, B, C, D, E, F)),
                    )

class Matrix_Chain_7_2_7(object):
    def __init__(self, n1 = 600, n2 = 200, n3 = 1500, n4 = 800, n5 = 1000):
        # equation: W = ABCDEF
        # A is lower triangular, 
        # D, E are upper triangular

        A = Matrix("A", (n1, n1), properties = [properties.INPUT, properties.LOWER_TRIANGULAR])
        B = Matrix("B", (n1, n2), properties = [properties.INPUT])
        C = Matrix("C", (n2, n3), properties = [properties.INPUT])
        D = Matrix("D", (n3, n3), properties = [properties.INPUT, properties.UPPER_TRIANGULAR])
        E = Matrix("E", (n3, n4), properties = [properties.INPUT, properties.UPPER_TRIANGULAR])
        F = Matrix("F", (n4, n5), properties = [properties.INPUT])

        W = Matrix("W", (n1, n5), properties = [properties.OUTPUT])

        self.eqns = Equations(
                            Equal(W, Times(
                                    Inverse(A),
                                    B, C,
                                    InverseTranspose(D),
                                    E, F
                            )),
                    )

class Properties_7_2_8(object):
    def __init__(self, n = 1000):
        # equation: X = (L1^-1 L2 U1^T + L3)^-1 U2^T
        # LX is lower triangular, UX is upper triangular

        L1 = Matrix("L1", (n, n), properties = [properties.INPUT, properties.LOWER_TRIANGULAR])
        L2 = Matrix("L2", (n, n), properties = [properties.INPUT, properties.LOWER_TRIANGULAR])
        L3 = Matrix("L3", (n, n), properties = [properties.INPUT, properties.LOWER_TRIANGULAR])
        U1 = Matrix("U1", (n, n), properties = [properties.INPUT, properties.UPPER_TRIANGULAR])
        U2 = Matrix("U2", (n, n), properties = [properties.INPUT, properties.UPPER_TRIANGULAR])
        X = Matrix("X", (n, n), properties = [properties.OUTPUT])

        self.eqns = Equations(
                            Equal(X,
                                Times(
                                    Inverse(
                                        Plus(
                                            Times(
                                                Inverse(L1),
                                                L2,
                                                Transpose(U1)
                                            ),
                                            L3
                                        )
                                    ),
                                    Transpose(U2)
                                )
                            )   
                    )

class Transposed_Kernel_7_2_9(object):
    def __init__(self, m = 1000, n = 1000):
        # equation: X = x^T A
        # LX is lower triangular, UX is upper triangular

        A = Matrix("A", (m, n), properties = [properties.INPUT])
        x = Vector("x", (m, sONE), properties = [properties.INPUT])
        y = Vector("y", (sONE, n), properties = [properties.OUTPUT])

        self.eqns = Equations(
                            Equal(y,
                                Times(
                                    Transpose(x),
                                    A
                                )
                            )   
                    )

class Transposed_Kernel_7_2_10(object):
    def __init__(self, m = 100, n = 1000):
        # m << n
        # equation: X = ASB^T
        # A is nxn
        # S is nxn, symmetric
        # B is mxn

        A = Matrix("A", (n, n), properties = [properties.INPUT])
        S = Matrix("S", (n, n), properties = [properties.INPUT, properties.SYMMETRIC])
        B = Matrix("B", (m, n), properties = [properties.INPUT])
        X = Matrix("X", (n, m), properties = [properties.OUTPUT])

        self.eqns = Equations(
                            Equal(X,
                                Times(
                                    A,
                                    S,
                                    Transpose(B)
                                )
                            )   
                    )

class Simplification_7_2_11(object):
    def __init__(self, m = 1000, n = 1000):
        # equation: X = 1/3A + 2/3A + B
        # A is mxn
        # B is mxn
        # X is mxn

        A = Matrix("A", (m, n), properties = [properties.INPUT])
        B = Matrix("B", (m, n), properties = [properties.INPUT])
        X = Matrix("X", (m, n), properties = [properties.OUTPUT])

        self.eqns = Equations(
                            Equal(X,
                                Plus(
                                    Times(A, ConstantScalar(1.5)),
                                    Times(A, ConstantScalar(-0.5)),
                                    B
                                )
                            )   
                    )

class Simplification_7_2_12(object):
    def __init__(self, m = 1000, n = 1000):
        # equation: X = 1/3A + 2/3A + B
        # A is mxn
        # B is mxn
        # X is mxn

        A = Matrix("A", (m, n), properties = [properties.INPUT])
        B = Matrix("B", (m, n), properties = [properties.INPUT])
        X = Matrix("X", (m, n), properties = [properties.OUTPUT])

        self.eqns = Equations(
                            Equal(X,
                                Plus(
                                    Times(A, ConstantScalar(0.34)),
                                    Times(A, ConstantScalar(0.66)),
                                    B
                                )
                            )   
                    )


class Example121(object):
    def __init__(self, n1 = 10):

        # TAGS
        # conjugate prior


        S = Matrix("S", (n1, n1), properties = [properties.SPD])
        P = Matrix("P", (n1, n1), properties = [properties.SPD])

        x = Vector("x", (n1, 1))
        m = Vector("m", (n1, 1))
        y = Vector("y", (n1, 1))

        alpha = Scalar("alpha")

        self.eqns = Equations(
                            Equal(y, 
                                Times(
                                    Inverse(Plus(Inverse(S), Times(alpha, Inverse(P)))),
                                    Plus(Times(Inverse(S), m), Times(alpha, Inverse(P), x))
                                    )
                                ) 
                            )


class Example122(object):
    def __init__(self):
        pass

        # albataineh2014

        # L = 5
        # K = 30 to 50
        # G = 63 or 64
        # d = {0, 1, 2, 3, 4}

        # H: (G+L-1) x G, full rank

        # S: G x K block diagonal

        # b: K-d vector

        # n: (G+L-1)-d vector

        # C: GxG diagonal, CC^H = I, digonal elements are +/-1+/-j

        # missing: properties of Q

        # b = S^H H^H (sigma^2 H H^H + Q)^-1 r


class Example123(object):
    def __init__(self):

        # TAGS
        # noschese2016

        n = 1500 # 10
        m = 1000 # 5

        A = Matrix("A", (n, m))
        A.set_property(properties.FULL_RANK)
        A.set_property(properties.INPUT)

        b = Vector("b", (n, sONE))
        b.set_property(properties.INPUT)

        mu = Scalar("mu")

        I = IdentityMatrix(m, m)

        x = Vector("x", (m, sONE))
        x.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(x, 
                                Times(Inverse(Plus(Times(Transpose(A), A), Times(mu, I))), Transpose(A), b)
                                ) 
                            )



