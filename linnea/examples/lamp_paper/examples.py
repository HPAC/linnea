from ...algebra.expression import Symbol, Scalar, Vector, Matrix, ConstantScalar, \
                                Equal, Plus, Times, Transpose, Inverse, \
                                InverseTranspose, InverseConjugate, \
                                InverseConjugateTranspose, \
                                ConjugateTranspose, Index, IdentityMatrix, Identity

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
                                            Transpose(H),
                                            # ConjugateTranspose(H),
                                            H
                                        ),
                                        Times(
                                            sigma,
                                            I
                                        )
                                    )
                                ),
                                Transpose(H),
                                # ConjugateTranspose(H),
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
        y = Vector("y", (n, sONE), properties = [properties.INPUT])

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

class Rank_1_Tensor_Update_7_1_8(object):
    def __init__(self, N = 1000, msd = 1000, nsd = 1000):
        pass

class EnsembleKalmarFilter_7_1_9_1(object):
    def __init__(self, n = 1000, m = 1000):
        # rao2015
        # the first variant
        # Ki modeled as two input and output vars
        # equation: Ki_O = P_b * H^T * (H * P_b * H^T + R)^-1
        # x_a = x_b + Ki_I * (y - H_x_b)
        # P_a = (I - K*H)*P_b
        # TODO: R is symmetric positive semi-definite, which property?

        minusone = ConstantScalar(-1.0)

        Ki_O = Matrix("Ki_O", (n, m), properties = [properties.OUTPUT])
        Ki_I = Matrix("Ki_I", (n, m), properties = [properties.INPUT])
        P_b = Matrix("P_b", (n, n), properties = [properties.INPUT, properties.SPD])
        P_a = Matrix("P_a", (n, n), properties = [properties.OUTPUT])
        H = Matrix("H", (m, n), properties = [properties.INPUT])
        R = Matrix("R", (m, m), properties = [properties.INPUT, properties.SYMMETRIC])
        I = IdentityMatrix(n, n)

        x_a = Vector("x_a", (n, sONE), properties = [properties.OUTPUT])
        x_b = Vector("x_b", (n, sONE), properties = [properties.INPUT])
        y = Vector("y", (m, sONE), properties = [properties.INPUT])
        H_x_b = Vector("H_x_b", (m, sONE), properties = [properties.INPUT])

        derivation.special_properties.add_expression(Plus(Times(H, P_b, Transpose(H)), R), {properties.SPD})

        self.eqns = Equations(
                            Equal(
                                Ki_O,
                                Times(
                                    P_b,
                                    Transpose(H),
                                    Inverse(
                                        Plus(
                                            Times(
                                                H,
                                                P_b,
                                                Transpose(H)
                                            ),
                                            R
                                        )
                                    )
                                )
                            ),
                            Equal(
                                x_a,
                                Plus(
                                    x_b,
                                    Times(
                                        Ki_I,
                                        Plus(
                                            y,
                                            Times(
                                                minusone,
                                                H_x_b
                                            )
                                        )
                                    )
                                )
                            ),
                            Equal(
                                P_a,
                                Times(
                                    Plus(
                                        I,
                                        Times(
                                            minusone,
                                            Ki_I,
                                            H
                                        )
                                    ),
                                    P_b
                                )
                            )
                    )


class EnsembleKalmarFilter_7_1_9_2(object):
    def __init__(self, n = 1000, m = 1000):
        # rao2015
        # the first variant
        # Ki modeled as two input and output vars
        # equation: Si_O = ((n-1)*I + Y^T * R^-1 * Y)^-1
        # K = X*Si_I*Y^T * R^-1
        # x_a = x_b + X*S*Y^T * R^-1 * (y - H_x_b)
        # P_a = X * S * X^T
        # TODO: R is symmetric positive semi-definite, which property?

        minusone = ConstantScalar(-1.0)
        n_scalar = ConstantScalar(n)

        K = Matrix("K", (n, m), properties = [properties.OUTPUT])
        Si_O = Matrix("Si_O", (n, m), properties = [properties.OUTPUT])
        Si_I = Matrix("Si_I", (n, m), properties = [properties.INPUT])
        X = Matrix("X", (m, n), properties = [properties.INPUT, properties.SPD])
        Y = Matrix("Y", (n, n), properties = [properties.INPUT, properties.SPD])
        R = Matrix("R", (m, m), properties = [properties.INPUT, properties.SYMMETRIC])
        I = IdentityMatrix(n, n)
        P_a = Matrix("P_a", (n, n), properties = [properties.OUTPUT])

        x_a = Vector("x_a", (n, sONE), properties = [properties.OUTPUT])
        x_b = Vector("x_b", (n, sONE), properties = [properties.INPUT])
        y = Vector("y", (m, sONE), properties = [properties.INPUT])
        H_x_b = Vector("H_x_b", (m, sONE), properties = [properties.INPUT])

        derivation.special_properties.add_expression(Plus(Times(Plus(n_scalar, minusone), I), Times(Transpose(Y), Inverse(R), Y)), {properties.SPD})

        self.eqns = Equations(
                            # Si_O = ((n-1)*I + Y^T * R^-1 * Y)^-1
                            Equal(
                                Si_O,
                                Inverse( Plus(
                                    Times(
                                        Plus(
                                            n_scalar,
                                            minusone
                                        ),
                                        I
                                    ),
                                    Times(
                                        Transpose(Y),
                                        Inverse(R),
                                        Y
                                    )
                                ))
                            ),
                            # K = X*Si_I*Y^T * R^-1
                            Equal(
                                K,
                                Times(
                                    X,
                                    Si_I,
                                    Transpose(Y),
                                    Inverse(R)
                                )
                            ),
                            # x_a = x_b + X*S*Y^T * R^-1 * (y - H_x_b)
                            Equal(
                                x_a,
                                Plus(
                                    x_b,
                                    Times(
                                        X,
                                        Si_I,
                                        Transpose(Y),
                                        Inverse(R),
                                        Plus(
                                            y,
                                            Times(
                                                minusone,
                                                H_x_b
                                            )
                                        )
                                    )
                                )
                            ),
                            # P_a = X * S * X^T
                            Equal(
                                P_a,
                                Times(
                                    X,
                                    Si_I,
                                    Transpose(X)
                                )
                            )
                    )

class SPA_7_1_12(object):
    def __init__(self, d = 1000, m = 1000, k = 8000, q = 1):
        # Trier2017
        # m < n
        # Algorithm nr 1 P^P
        # Y = (AA^T)^q A_I
        # B = QQ^T*A


        A = Matrix("A", (d, m), properties = [properties.INPUT])
        A_I = Matrix("A_I", (d, k), properties = [properties.INPUT])
        Q = Matrix("Q", (d, k), properties = [properties.INPUT, properties.ORTHOGONAL])
        B = Matrix("B", (d, m), properties = [properties.OUTPUT])
        Y = Matrix("Y", (d, k), properties = [properties.OUTPUT])

        power_arg = Times(A, Transpose(A))
        if q == 0:
            power_expr = Identity()
        else:
            power_expr = power_arg
            for i in range(1, q):
                power_expr = Times(power_expr, power_arg)

        self.eqns = Equations(
                            Equal(
                                Y,
                                Times(
                                    power_expr,
                                    A_I
                                )
                            ),
                            Equal(
                                B,
                                Times(
                                    Q,
                                    Transpose(Q),
                                    A
                                )
                            )
                    )

class ImageRestoration_7_1_13_1(object):
    def __init__(self, n = 1500, m = 1000):
        # Trier2017
        # m < n
        # Algorithm nr 1 P^P
        # x = (H^t * H + lambda * sigma^2 * I_n)^-1 * (H^T * y + lambda * sigma^2 * (v - u))

        minusone = ConstantScalar(-1)
        lambda_ = Scalar("lambda")
        sigma_ = Scalar("sigma_sq")

        H = Matrix("H", (m, n), properties = [properties.INPUT, properties.FULL_RANK])
        I = IdentityMatrix(n, n)

        v_k = Vector("v_k", (n, sONE), properties = [properties.INPUT])
        u_k = Vector("u_k", (n, sONE), properties = [properties.INPUT])
        y = Vector("y", (m, sONE), properties = [properties.INPUT])
        x = Vector("x", (n, sONE), properties = [properties.OUTPUT])


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
                                                    minusone,
                                                    u_k
                                                )
                                            )
                                        )
                                    )
                                )
                            )
                    )

class ImageRestoration_7_1_13_2(object):
    def __init__(self, n = 1500, m = 1000, single = False):
        # Trier2017
        # m < n
        # Algorithm nr 2
        # H^dag = H^T (HH^T)^-1
        # y_k = H^dag y + (I - H^dag H)x_k
        # H_dag is modeled input/output

        minusone = ConstantScalar(-1.0)
        lambda_ = Scalar("lambda")
        sigma_ = Scalar("sigma_sq")

        H = Matrix("H", (m, n), properties = [properties.INPUT, properties.FULL_RANK])
        H_dag_I = Matrix("H_dag_I", (n, m), properties = [properties.INPUT])
        H_dag_O = Matrix("H_dag_O", (n, m), properties = [properties.OUTPUT])
        I = IdentityMatrix(n, n)

        y_k = Vector("y_k", (n, sONE), properties = [properties.OUTPUT])
        y = Vector("y", (m, sONE), properties = [properties.INPUT])
        x = Vector("x", (n, sONE), properties = [properties.INPUT])

        h_dag = Times(
                    Transpose(H),
                    Inverse(
                        Times(
                            H,
                            Transpose(H)
                        )
                    )
                )
        if single:
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
                                                minusone,
                                                h_dag,
                                                H
                                            )
                                        ),
                                        x
                                    )
                                )
                            )
                        )
        else:
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
                                                minusone,
                                                H_dag_I,
                                                H
                                            )
                                        ),
                                        x
                                    )
                                )
                            )
                        )

class Tikhonov_7_1_14(object):
    def __init__(self, n = 1500, m = 1000):

        # TAGS
        # noschese2016
        # former Example123
        # eq: x = (A^T A + \mu^2 I)^-1 A^T b
        # eq. 1.5, 2.7

        A = Matrix("A", (n, m), properties = [properties.INPUT, properties.FULL_RANK])
        I = IdentityMatrix(m, m)
        b = Vector("b", (n, sONE), properties = [properties.INPUT])
        x = Vector("x", (m, sONE), properties = [properties.OUTPUT])
        mu = Scalar("mu_sq")

        derivation.special_properties.add_expression(Plus(Times(Transpose(A), A), Times(mu, I)), {properties.SPD})

        self.eqns = Equations(
                            Equal(x,
                                Times(Inverse(Plus(Times(Transpose(A), A), Times(mu, I))), Transpose(A), b)
                                )
                            )

class CDMA_7_1_15(object):
    def __init__(self, G=63, K=30, L=5):
        pass

        # albataineh2014

        # L = 5
        # K = 30 to 50
        # G = 63 or 64
        # d = {0, 1, 2, 3, 4}
        # G1 >= G + L -1
        # H: G1 x G, full rank

        # S: G x K block diagonal

        # b: K-d vector

        # n: G1-d vector

        # C: GxG diagonal, CC^H = I, digonal elements are +/-1+/-j

        # missing: properties of Q
        # "and n is the (G1) dimensional channel noise vector with covariance matrix, say, Q"
        # covariances matrices are symmetric and positive semi-definite 

        # b = S^H H^H (sigma^2 H H^H + Q)^-1 r
        G1 = G + L -1

        sigma_sq = Scalar("sigma_sq")
        S = Matrix("S", (G, K), properties = [properties.INPUT])
        H = Matrix("H", (G1, G), properties = [properties.INPUT, properties.FULL_RANK])
        Q = Matrix("Q", (G1, G1), properties = [properties.INPUT, properties.SYMMETRIC])
        r = Vector("r", (G1, sONE), properties = [properties.INPUT])
        b = Vector("b", (K, sONE), properties = [properties.OUTPUT])

        self.eqns = Equations(
                            Equal(b,
                                Times(
                                    Transpose(S),
                                    Transpose(H),                                    
                                    # ConjugateTranspose(S),
                                    # ConjugateTranspose(H),
                                    Inverse(Plus(
                                        Times(
                                            sigma_sq,
                                            H,
                                            Transpose(H)
                                            # ConjugateTranspose(H)
                                        ),
                                        Q
                                    )),
                                    r
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





