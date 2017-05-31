

from .algebra.expression import Symbol, Scalar, Vector, Matrix, ConstantScalar, \
                                Equal, Plus, Times, Transpose, Inverse, \
                                InverseTranspose, InverseConjugate, \
                                InverseConjugateTranspose, \
                                ConjugateTranspose, Index, IdentityMatrix

from .algebra.properties import Property as properties

from .algebra.equations import Equations

from . import derivation

sONE = 1
sZERO = 0

#######################
### Diego's examples

class Example001(object):
    def __init__(self):

        # TAGS
        # matrix chain

        n = 10

        x = Vector("x", (n, sONE))
        x.set_property(properties.INPUT)

        y = Vector("y", (n, sONE))
        y.set_property(properties.INPUT)

        z = Vector("z", (n, sONE))
        z.set_property(properties.INPUT)

        alpha = Scalar("alpha")
        alpha.set_property(properties.OUTPUT)

        # alpha = x^T z x^T y
        self.eqns = Equations(Equal(alpha, Times(Transpose(x), z, Transpose(x), y)))

        

###

class Example002(object):
    def __init__(self):

        # TAGS
        # SOP, POS, CSE

        n = 10

        X = Matrix("X", (n, n))
        X.set_property(properties.INPUT)

        Y = Matrix("Y", (n, n))
        Y.set_property(properties.INPUT)

        v = Vector("v", (n, sONE))
        v.set_property(properties.INPUT)

        w = Vector("w", (n, sONE))
        w.set_property(properties.OUTPUT)

        # w = X Y^T v + Y X^T v
        self.eqns = Equations(Equal(w, Plus(Times(X, Transpose(Y), v), Times(Y, Transpose(X), v))))

###

class Example003(object):
    def __init__(self):

        # TAGS
        # matrix chain

        n = 10

        L = Matrix("L", (n, n))
        L.set_property(properties.INPUT)
        L.set_property(properties.LOWER_TRIANGULAR)

        v = Vector("v", (n, sONE))
        v.set_property(properties.INPUT)

        u = Vector("u", (n, sONE))
        u.set_property(properties.INPUT)

        alpha = Scalar("alpha")
        alpha.set_property(properties.OUTPUT)

        # alpha = u^T L^-1 L^-1 v
        self.eqns = Equations(Equal(alpha, Times(Transpose(u), Inverse(L), Inverse(L), v)))

        expr = Times(Transpose(u), Inverse(L), Inverse(L), v)

###

class Example004(object):
    def __init__(self):

        # TAGS
        # matrix chain, CSE

        n = 10

        x = Vector("x", (n, sONE))
        x.set_property(properties.INPUT)

        y = Vector("y", (n, sONE))
        y.set_property(properties.INPUT)

        alpha = Scalar("alpha")
        alpha.set_property(properties.OUTPUT)

        # alpha = x^T y x^T y
        self.eqns = Equations(Equal(alpha, Times(Transpose(x), y, Transpose(x), y)))

###

class Example005(object):
    def __init__(self):

        # TAGS
        # solve equation

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.INPUT)
        A.set_property(properties.SPD)
        A.set_property(properties.SQUARE)
        A.set_property(properties.NON_SINGULAR)

        B = Matrix("B", (n, n))
        B.set_property(properties.INPUT)

        X = Matrix("X4", (n, n))
        X.set_property(properties.OUTPUT)

        # A X = B
        self.eqns = Equations(Equal(Times(A, X), B))

###

class Example006(object):
    def __init__(self):

        # TAGS
        # explicit inverse, matrix chain

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.INPUT)
        A.set_property(properties.SYMMETRIC)
        A.set_property(properties.SQUARE)
        A.set_property(properties.NON_SINGULAR)

        B = Matrix("B", (n, n))
        B.set_property(properties.OUTPUT)
        B.set_property(properties.SYMMETRIC)
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)

        # B = A^-1
        self.eqns = Equations(Equal(B, Inverse(A)))

###

class Example007(object):
    def __init__(self):

        # TAGS
        # least squares, QR

        n = 1500
        m = 1000

        X = Matrix("X", (n, m))
        X.set_property(properties.INPUT)
        X.set_property(properties.FULL_RANK)

        y = Vector("y", (n, sONE))
        y.set_property(properties.INPUT)

        b = Vector("b", (m, sONE))
        b.set_property(properties.OUTPUT)

        # b = (X^T X)^-1 X^T y
        self.eqns = Equations(Equal(b, Times(Inverse(Times(Transpose(X), X)), Transpose(X), y)))

###

class Example008(object):
    def __init__(self):

        # TAGS
        # simple

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.INPUT)

        B = Matrix("B", (n, n))
        B.set_property(properties.INPUT)

        C = Matrix("C", (n, n))
        C.set_property(properties.INPUT)

        D = Matrix("D", (n, n))
        D.set_property(properties.INPUT)

        E = Matrix("E", (n, n))
        E.set_property(properties.INPUT)

        X = Matrix("X", (n, n))
        X.set_property(properties.OUTPUT)

        # X = A B + C D + E
        self.eqns = Equations(Equal(X, Plus(Times(A, B), Times(C, D), E)))

###

class Example009(object):
    def __init__(self):

        # TAGS
        # CSE, matrix chain

        n = 10

        L = Matrix("L", (n, n))
        L.set_property(properties.INPUT)
        L.set_property(properties.LOWER_TRIANGULAR)

        A = Matrix("A", (n, n))
        A.set_property(properties.INPUT)

        X = Matrix("X", (n, n))
        X.set_property(properties.OUTPUT)

        # X = A^T L^-1 L^-T A
        self.eqns = Equations(Equal(X, Times(Transpose(A), Inverse(L), Transpose(Inverse(L)), A)))

### GWAS in Demo-Rev missing and derivative

###############################
### my examples


class Example010(object):
    def __init__(self):

        # TAGS
        # complex inverse

        n = 10
        m = 5

        S = Matrix("S", (n, n))
        S.set_property(properties.SPD)
        S.set_property(properties.INPUT)

        CP = Matrix("CP", (n, m))
        CP.set_property(properties.COLUMN_PANEL)
        CP.set_property(properties.FULL_RANK)
        CP.set_property(properties.INPUT)

        Y = Matrix("Y", (m, m))
        Y.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(Y, Inverse(Times(Transpose(CP), S, CP ) )))


class Example011(object):
    def __init__(self):

        # TAGS
        # symmetric product

        n = 10
        m = 5

        S = Matrix("S", (n, n))
        S.set_property(properties.SPD)
        S.set_property(properties.INPUT)

        CP = Matrix("CP", (n, m))
        CP.set_property(properties.COLUMN_PANEL)
        CP.set_property(properties.FULL_RANK)
        CP.set_property(properties.INPUT)

        Y = Matrix("Y", (m, m))
        Y.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(Y, Times(Transpose(CP), S, CP )))


class Example012(object):
    def __init__(self):
        
        # TAGS
        # simple

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        W.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(W, Times(Transpose(A), A)))


class Example013(object):
    def __init__(self):

        # TAGS
        # simple, matrix chain

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.INPUT)

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)
        B.set_property(properties.INPUT)

        C = Matrix("C", (n, n))
        C.set_property(properties.INPUT)

        D = Matrix("D", (n, n))
        D.set_property(properties.DIAGONAL)
        D.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        W.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(W, Times(A, B, C, D)))

class Example014(object):
    def __init__(self):

        # TAGS
        # matrix chain, CSE, symmetric product

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.INPUT)

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)
        B.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        W.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(W, Times(Transpose(B), Transpose(A), A, B)))


class Example015(object):
    def __init__(self):

        # TAGS
        # simple

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.INPUT)

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)
        B.set_property(properties.INPUT)

        C = Matrix("C", (n, n))
        C.set_property(properties.INPUT)

        D = Matrix("D", (n, n))
        D.set_property(properties.DIAGONAL)
        D.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        W.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(W, Plus(A, B, C, D)))


class Example016(object):
    def __init__(self):

        # TAGS
        # complex inverse, random

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.INPUT)

        D = Matrix("D", (n, n))
        D.set_property(properties.DIAGONAL)
        D.set_property(properties.INPUT)

        S = Matrix("S", (n, n))
        S.set_property(properties.SPD)
        S.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        W.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(W, Times(A, Inverse(S), Inverse(Plus(A, Inverse(S) ) ), D)))


class Example017(object):
    def __init__(self):

        # TAGS
        # complex inverse, random

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.INPUT)

        D = Matrix("D", (n, n))
        D.set_property(properties.DIAGONAL)
        D.set_property(properties.INPUT)

        S = Matrix("S", (n, n))
        S.set_property(properties.SPD)
        S.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        W.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(W, Times(A, Inverse(S), Inverse(Plus(A, Inverse(D) ) ), D)))


class Example018(object):
    def __init__(self):

        # TAGS
        # random, explicit inverse, sum

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.INPUT)

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)
        B.set_property(properties.INPUT)

        C = Matrix("C", (n, n))
        C.set_property(properties.INPUT)

        D = Matrix("D", (n, n))
        D.set_property(properties.DIAGONAL)
        D.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        W.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(W, Plus(A, Inverse(B), Inverse(C), D)))


class Example019(object):
    def __init__(self):

        # TAGS
        # random, explicit inverse, sum

        n = 10

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)
        B.set_property(properties.INPUT)

        C = Matrix("C", (n, n))
        C.set_property(properties.INPUT)

        D = Matrix("D", (n, n))
        D.set_property(properties.DIAGONAL)
        D.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        W.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(W, Plus(Inverse(B), Inverse(C), D)))


class Example020(object):
    def __init__(self):

        # TAGS
        # explicit inverse

        n = 10

        S = Matrix("S", (n, n))
        S.set_property(properties.SPD)
        S.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        W.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(W, Inverse(S)))


class Example021(object):
    def __init__(self):

        # TAGS
        # matrix chain, linear system

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.INPUT)

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)
        # B.set_property(properties.SPD)
        B.set_property(properties.INPUT)

        x = Vector("x", (n, sONE))
        x.set_property(properties.INPUT)

        w = Vector("w", (n, sONE))
        w.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(w, Times(A, Inverse(B), x)))


class Example022(object):
    def __init__(self):

        # TAGS
        # linear system, matrix chain

        n = 10

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)
        B.set_property(properties.INPUT)

        x = Vector("x", (n, sONE))
        x.set_property(properties.INPUT)

        w = Vector("w", (n, sONE))
        w.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(w, Times(Inverse(B), x)))


class Example023(object):
    def __init__(self):

        # TAGS
        # triangular system

        n = 10

        L = Matrix("L", (n, n))
        L.set_property(properties.LOWER_TRIANGULAR)
        L.set_property(properties.NON_SINGULAR)
        L.set_property(properties.INPUT)

        x = Vector("x", (n, sONE))
        x.set_property(properties.INPUT)

        w = Vector("w", (n, sONE))
        w.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(w, Times(Inverse(L), x)))


class Example024(object):
    def __init__(self):

        # TAGS
        # SPD system

        n = 10

        S = Matrix("S", (n, n))
        S.set_property(properties.SPD)
        S.set_property(properties.INPUT)

        x = Vector("x", (n, sONE))
        x.set_property(properties.INPUT)

        w = Vector("w", (n, sONE))
        w.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(w, Times(Inverse(S), x)))


class Example025(object):
    def __init__(self):

        # TAGS
        # CSE, SOP, POS, random

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.INPUT)

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)
        B.set_property(properties.INPUT)

        S = Matrix("S", (n, n))
        S.set_property(properties.SPD)
        S.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        W.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(W, Plus(Times(A, B, S), Times(Transpose(B), Transpose(A), S) )))


class Example026(object):
    def __init__(self):

        # TAGS
        # random, sum, simple

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.INPUT)

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)
        B.set_property(properties.INPUT)

        C = Matrix("C", (n, n))
        C.set_property(properties.INPUT)

        D = Matrix("D", (n, n))
        D.set_property(properties.DIAGONAL)
        D.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        W.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(W, Plus(A, Transpose(B), C, D)))


class Example027(object):
    def __init__(self):

        # TAGS
        # transposed kernel

        n = 10

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)
        B.set_property(properties.INPUT)

        C = Matrix("C", (n, n))
        C.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        W.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(W, Times(Transpose(B), Transpose(C))))


class Example028(object):
    def __init__(self):

        # TAGS
        # transposed kernel, CSE

        n = 10

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)
        B.set_property(properties.INPUT)

        C = Matrix("C", (n, n))
        C.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        W.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(W, Times(Transpose(B), Transpose(C), C, B)))


class Example029(object):
    def __init__(self):

        # TAGS
        # transposed kernel

        n = 10

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)
        B.set_property(properties.INPUT)

        x = Vector("x", (n, sONE))
        x.set_property(properties.INPUT)

        y = Vector("y", (n, sONE))
        y.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        W.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(W, Times(y, Transpose(x), B)))


class Example030(object):
    def __init__(self):

        # TAGS
        # transposed kernel

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.INPUT)

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)
        B.set_property(properties.INPUT)

        x = Vector("x", (n, sONE))
        x.set_property(properties.INPUT)

        y = Vector("y", (n, sONE))
        y.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        W.set_property(properties.OUTPUT)

        U = Matrix("U", (n, n))
        U.set_property(properties.OUTPUT)

        self.eqns = Equations(
                        Equal(W, Times(y, Transpose(x), B)),
                        Equal(U, Times(Transpose(A), A))
                        )

class Example031(object):
    def __init__(self):

        # TAGS
        # blocked algorithm

        n = 10
        m = 10
        k = 10

        L00 = Matrix("L00", (n, n))
        L00.set_property(properties.INPUT)
        L00.set_property(properties.LOWER_TRIANGULAR)

        L11 = Matrix("L11", (m, m))
        L11.set_property(properties.INPUT)
        L11.set_property(properties.LOWER_TRIANGULAR)

        L22 = Matrix("L22", (k, k))
        L22.set_property(properties.INPUT)
        L22.set_property(properties.LOWER_TRIANGULAR)

        L10 = Matrix("L10", (m, n))
        L10.set_property(properties.INPUT)

        L21 = Matrix("L21", (k, m))
        L21.set_property(properties.INPUT)

        X10_in = Matrix("X10_in", (m, n))
        X10_in.set_property(properties.INPUT)

        X10 = Matrix("X10", (m, n))
        X10.set_property(properties.OUTPUT)

        X11 = Matrix("X11", (m, m))
        X11.set_property(properties.OUTPUT)

        X20_in = Matrix("X20_in", (k, n))
        X20_in.set_property(properties.INPUT)

        X20 = Matrix("X20", (k, n))
        X20.set_property(properties.OUTPUT)

        X21 = Matrix("X21", (k, m))
        X21.set_property(properties.OUTPUT)

        self.eqns = Equations(
                        Equal(X10, Times(X10_in, Inverse(L00))),
                        Equal(X20, Plus(X20_in, Times(Inverse(L22), L21, Inverse(L11), L10))),
                        Equal(X11, Inverse(L11)),
                        # Equal(X21, Minus([Times(Inverse(L22), L21))), # this is the correct equation, but minus is a problem
                        Equal(X21, Times(Inverse(L22), L21)),
                        )
        

class Example032(object):
    def __init__(self):

        # TAGS
        # CSE times

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.INPUT)

        B = Matrix("B", (n, n))
        B.set_property(properties.INPUT)

        C = Matrix("C", (n, n))
        C.set_property(properties.INPUT)

        D = Matrix("D", (n, n))
        D.set_property(properties.INPUT)

        X = Matrix("X", (n, n))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n, n))
        Y.set_property(properties.OUTPUT)

        self.eqns = Equations(
                        Equal(X, Times(A, B, C)),
                        Equal(Y, Times(B, C, D))
                        )


class Example033(object):
    def __init__(self):

        # TAGS
        # CSE plus

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.INPUT)

        B = Matrix("B", (n, n))
        B.set_property(properties.INPUT)

        C = Matrix("C", (n, n))
        C.set_property(properties.INPUT)

        D = Matrix("D", (n, n))
        D.set_property(properties.INPUT)

        X = Matrix("X", (n, n))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n, n))
        Y.set_property(properties.OUTPUT)

        self.eqns = Equations(
                        Equal(X, Plus(A, B, C)),
                        Equal(Y, Plus(B, C, D))
                        )


class Example034(object):
    def __init__(self):

        # TAGS
        # CSE plus

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.INPUT)

        B = Matrix("B", (n, n))
        B.set_property(properties.INPUT)

        C = Matrix("C", (n, n))
        C.set_property(properties.INPUT)

        D = Matrix("D", (n, n))
        D.set_property(properties.INPUT)

        X = Matrix("X", (n, n))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n, n))
        Y.set_property(properties.OUTPUT)

        Z = Matrix("Z", (n, n))
        Z.set_property(properties.AUXILIARY)

        self.eqns = Equations(
                        Equal(Z, Plus(A, B)),
                        Equal(X, Plus(C, Z)),
                        Equal(Y, Plus(Z, D))
                        )


class Example035(object):
    def __init__(self):

        # TAGS
        # matrix chain, simple

        n1 = 10
        n2 = 20
        n3 = 15
        n4 = 30

        A = Matrix("A", (n1, n2))
        A.set_property(properties.INPUT)

        B = Matrix("B", (n2, n3))
        B.set_property(properties.INPUT)

        C = Matrix("C", (n3, n4))
        C.set_property(properties.INPUT)

        d = Vector("d", (n4, sONE))
        d.set_property(properties.INPUT)

        x = Vector("x", (n1, sONE))
        x.set_property(properties.OUTPUT)

        self.eqns = Equations(
                        Equal(x, Times(A, B, C, d))
                        )
        

class Example036(object):
    def __init__(self):

        # TAGS
        # random, complex inverse

        n1 = 10
        n2 = 20
        n3 = 15

        A = Matrix("A", (n1, n1))
        A.set_property(properties.INPUT)
        A.set_property(properties.SPD)

        B = Matrix("B", (n1, n2))
        B.set_property(properties.INPUT)

        C = Matrix("C", (n2, n1))
        C.set_property(properties.INPUT)

        D = Matrix("D", (n1, n3))
        D.set_property(properties.INPUT)

        X = Matrix("X", (n1, n3))
        X.set_property(properties.OUTPUT)

        self.eqns = Equations(
                        Equal(X, Times(Inverse(Times(A, B, C)), D))
                        )


class Example037(object):
    def __init__(self):

        # TAGS
        # random, simple

        n1 = 10
        n2 = 20
        n3 = 15
        n4 = 30
        n5 = 10
        n6 = 20

        A = Matrix("A", (n1, n2))
        A.set_property(properties.INPUT)

        B = Matrix("B", (n2, n3))
        B.set_property(properties.INPUT)

        C = Matrix("C", (n3, n4))
        C.set_property(properties.INPUT)

        d = Vector("d", (n4, sONE))
        d.set_property(properties.INPUT)

        x = Vector("x", (n1, sONE))
        x.set_property(properties.OUTPUT)

        E = Matrix("E", (n5, n6))
        E.set_property(properties.INPUT)

        F = Matrix("F", (n5, n6))
        F.set_property(properties.INPUT)

        G = Matrix("G", (n5, n6))
        G.set_property(properties.INPUT)

        Y = Matrix("Y", (n5, n6))
        Y.set_property(properties.OUTPUT)

        self.eqns = Equations(
                        Equal(x, Times(A, B, C, d)),
                        Equal(Y, Plus(E, F, G))
                        )


class Example038(object):
    def __init__(self):

        # TAGS
        # indexed matrix chain, indices

        i = Index("i", 2)
        j = Index("j", 20)
        k = Index("k", 100)
        n1 = 10
        n2 = 10
        n3 = 10
        n4 = 10
        n5 = 10
        n6 = 10
        n7 = 10

        A = Matrix("A", (n1, n2), set())
        B = Matrix("B", (n2, n3), set())
        C = Matrix("C", (n3, n4), set([i]))
        D = Matrix("D", (n4, n5), set([i, j]))
        E = Matrix("E", (n5, n6), set([i, j, k]))
        F = Matrix("F", (n6, n7), set())
        # G = Matrix("G", (10, 10), set([]))
        # H = Matrix("H", (10, 10), set([]))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)
        C.set_property(properties.INPUT)
        D.set_property(properties.INPUT)
        E.set_property(properties.INPUT)
        F.set_property(properties.INPUT)
        # G.set_property(properties.INPUT)
        # H.set_property(properties.INPUT)

        X = Matrix("X", (n1, n7), set([i, j ,k]))

        X.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            # Equal(X, Times(A, B, C, D, E, F, G, H))
                            Equal(X, Times(A, B, C, D, E, F)) 
                            )




class Example039(object):
    def __init__(self):
        
        # TAGS
        # indexed matrix sum, indices, number of entries

        i = Index("i", 2)
        j = Index("j", 20)
        k = Index("k", 100)
        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n2), set([]))
        B = Matrix("B", (n1, n2), set([]))
        C = Matrix("C", (n1, n2), set([i]))
        D = Matrix("D", (n1, n2), set([i, j]))
        E = Matrix("E", (n1, n2), set([i, j ,k]))
        F = Matrix("F", (n1, n2), set([]))
        # G = Matrix("G", (n1, n2), set([]))
        # H = Matrix("H", (n1, n2), set([]))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)
        C.set_property(properties.INPUT)
        D.set_property(properties.INPUT)
        E.set_property(properties.INPUT)
        F.set_property(properties.INPUT)
        # G.set_property(properties.INPUT)
        # H.set_property(properties.INPUT)

        A.set_property(properties.LOWER_TRIANGULAR)
        B.set_property(properties.LOWER_TRIANGULAR)
        F.set_property(properties.LOWER_TRIANGULAR)


        X = Matrix("X", (n1, n2), set([i, j ,k]))
        X.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            # Equal(X, Plus(A, B, C, D, E, F, G, H)) 
                            Equal(X, Plus(A, B, C, D, E, F)) 
                            )
        

class Example040(object):
    def __init__(self):

        # TAGS
        # CSE times

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)
        C.set_property(properties.INPUT)

        X = Matrix("X", (n1, n1))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n1, n1))
        Y.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(X, Plus(A, B, C, Transpose(A), Transpose(B))),
                            Equal(Y, Plus(A, B))
                            )


class Example041(object):
    def __init__(self):

        # TAGS
        # CSE plus

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        D = Matrix("D", (n1, n1))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)
        C.set_property(properties.INPUT)
        D.set_property(properties.INPUT)

        X = Matrix("X", (n1, n1))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n1, n1))
        Y.set_property(properties.OUTPUT)

        Z = Matrix("Z", (n1, n1))
        Z.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(X, Plus(Transpose(A), B, C,)),
                            Equal(Y, Plus(B, C, D)),
                            Equal(Z, Plus(Transpose(B), A, C))
                            )


class Example042(object):
    def __init__(self):

        # TAGS
        # GWAS, generalized least squares (inverse missing), complex inverse
        # symmetric product

        n = 10
        m = 5

        S = Matrix("S", (n, n))
        S.set_property(properties.SPD)
        S.set_property(properties.INPUT)

        CP = Matrix("CP", (n, m))
        CP.set_property(properties.COLUMN_PANEL)
        CP.set_property(properties.FULL_RANK)
        CP.set_property(properties.INPUT)

        z = Vector("z", (m, sONE))
        z.set_property(properties.OUTPUT)

        y = Vector("y", (n, sONE))
        y.set_property(properties.INPUT)

        self.eqns = Equations(Equal(z, Times(Inverse(Times(Transpose(CP), S, CP ) ), Transpose(CP), Inverse(S), y)))


class Example043(object):
    def __init__(self):

        # TAGS
        # CSE times

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))

        A.set_property(properties.INPUT)
        A.set_property(properties.LOWER_TRIANGULAR)
        A.set_property(properties.SQUARE)
        B.set_property(properties.INPUT)
        B.set_property(properties.LOWER_TRIANGULAR)
        B.set_property(properties.SQUARE)
        C.set_property(properties.INPUT)

        X = Matrix("X", (n1, n1))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n1, n1))
        Y.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(X, Times(A, B, C, Transpose(Inverse(A)), Transpose(Inverse(B)))),
                            Equal(Y, Times(A, B))
                            )


class Example044(object):
    def __init__(self):

        # TAGS
        # CSE times

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))

        A.set_property(properties.INPUT)
        A.set_property(properties.LOWER_TRIANGULAR)
        A.set_property(properties.SQUARE)
        B.set_property(properties.INPUT)
        B.set_property(properties.LOWER_TRIANGULAR)
        B.set_property(properties.SQUARE)
        C.set_property(properties.INPUT)

        X = Matrix("X", (n1, n1))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n1, n1))
        Y.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(X, Times(Transpose(B), Transpose(A), C, Transpose(Inverse(A)), Transpose(Inverse(B)))),
                            Equal(Y, Times(A, B))
                            )

class Example045(object):
    def __init__(self):

        # TAGS
        # CSE times

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)
        C.set_property(properties.INPUT)

        X = Matrix("X", (n1, n1))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n1, n1))
        Y.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(X, Times(A, B, C, Transpose(Inverse(A)), Transpose(Inverse(B)))),
                            Equal(Y, Times(A, B))
                            )


class Example046(object):
    def __init__(self):

        # TAGS
        # CSE times

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)
        C.set_property(properties.INPUT)

        X = Matrix("X", (n1, n1))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n1, n1))
        Y.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(X, Times(A, B, C, Transpose(Inverse(A)), Transpose(Inverse(B)))),
                            Equal(Y, Times(A, B))
                            )


class Example047(object):
    def __init__(self):

        # TAGS
        # POS, SOP

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        D = Matrix("D", (n1, n1))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)
        C.set_property(properties.INPUT)
        D.set_property(properties.INPUT)

        X = Matrix("X", (n1, n1))
        X.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            # Equal(X, Times(Plus(A, B), Plus(C, Times(D, Plus(B, C)))))
                            # Equal(X, Times(Plus(A, B), Plus(C, D)))
                            # Equal(X, Times(Plus(A, B), C))
                            Equal(X, Plus(Times(A, Plus(C, D)), Times(B, Plus(C, D))))
                            )


class Example048(object):
    def __init__(self):

        # TAGS
        # POS, SOP

        n1 = 10
        n2 = 10 

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        D = Matrix("D", (n1, n1))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)
        C.set_property(properties.INPUT)
        D.set_property(properties.INPUT)

        X = Matrix("X", (n1, n1))
        X.set_property(properties.OUTPUT)

        # from patternmatcher.expression import to_SOP

        self.eqns = Equations(
                            # Equal(X, Times(Plus(A, B), Plus(C, Times(D, Plus(B, C)))))
                            Equal(X, Times(Plus(A, B), Plus(C, D)))
                            # Equal(X, to_SOP(Times(Plus(A, Transpose(Inverse(Plus(A, B)))), Plus(C, Inverse(Plus(A, B))))))
                            # Equal(X, Times(Plus(A, B), C))
                            # Equal(X, Plus(Times(A, Plus(B, C)), Times(D, Plus(B, C))))
                            )


class Example049(object):
    def __init__(self):

        # TAGS
        # linear system, matrix chain

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))

        A.set_property(properties.INPUT)
        A.set_property(properties.NON_SINGULAR)
        B.set_property(properties.INPUT)

        X = Matrix("X", (n1, n1))
        X.set_property(properties.OUTPUT)

        # from patternmatcher.expression import to_SOP

        self.eqns = Equations(
                            Equal(X, Times(InverseTranspose(A), B))
                            )


class Example050(object):
    def __init__(self):

        # TAGS
        # CSE plus times, random

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)
        C.set_property(properties.INPUT)

        X = Matrix("X", (n1, n1))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n1, n1))
        Y.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(X, Plus(Times(A, B), C, Inverse(Times(A, B))))
                            )


class Example051(object):
    def __init__(self):

        # TAGS
        # CSE, random, hard

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)
        C.set_property(properties.INPUT)

        X = Matrix("X", (n1, n1))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n1, n1))
        Y.set_property(properties.OUTPUT)

        Z = Matrix("Z", (n1, n1))
        Z.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(X, Plus(A, C, Inverse(Times(A, B)))),
                            Equal(Y, Plus(B, C, Inverse(Times(A, B)))),
                            Equal(Z, Plus(B, C, InverseTranspose(Times(A, B))))
                            )


class Example052(object):
    def __init__(self):

        # TAGS
        # CSE plus times, random, hard

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)
        C.set_property(properties.INPUT)

        X = Matrix("X", (n1, n1))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n1, n1))
        Y.set_property(properties.OUTPUT)

        Z = Matrix("Z", (n1, n1))
        Z.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(X, Plus(A, C, Inverse(Times(A, B)))),
                            Equal(Y, Plus(A, C, Inverse(Times(A, C)))),
                            Equal(Y, Plus(B, C, Inverse(Times(A, C)))),
                            Equal(Z, Plus(B, C, InverseTranspose(Times(A, B))))
                            )


class Example053(object):
    def __init__(self):

        # TAGS
        # CSE times

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        D = Matrix("D", (n1, n1))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)
        C.set_property(properties.INPUT)
        D.set_property(properties.INPUT)

        X = Matrix("X", (n1, n1))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n1, n1))
        Y.set_property(properties.OUTPUT)

        Z = Matrix("Z", (n1, n1))
        Z.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(X, Times(A, B, C, D)),
                            Equal(Y, Times(C, D, A, B))
                            )

class Example054(object):
    def __init__(self):

        # TAGS
        # SOP, POS

        n1 = 1000
        n2 = 1000

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        D = Matrix("D", (n1, n1))
        E = Matrix("E", (n1, n1))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)
        C.set_property(properties.INPUT)
        D.set_property(properties.INPUT)
        E.set_property(properties.INPUT)

        X = Matrix("X", (n1, n1))
        X.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            # Equal(X, Times(Plus(A, B), Plus(C, Times(D, Plus(B, C)))))
                            # Equal(X, Times(Plus(A, B), Plus(C, D)))
                            # Equal(X, Times(Plus(A, B), C))
                            Equal(X, Plus(Times(A, Plus(C, D, E)), Times(B, Plus(C, D))))
                            )


class Example055(object):
    def __init__(self):

        # TAGS
        # CSE times, random

        n1 = 10
        n2 = 10
        n3 = 10

        A = Matrix("A", (n1, n1))
        A.set_property(properties.INPUT)
        A.set_property(properties.SPD)

        B = Matrix("B", (n1, n2))
        B.set_property(properties.INPUT)

        C = Matrix("C", (n2, n1))
        C.set_property(properties.INPUT)

        D = Matrix("D", (n1, n3))
        D.set_property(properties.INPUT)

        X = Matrix("X", (n1, n3))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n1, n1))
        Y.set_property(properties.OUTPUT)

        self.eqns = Equations(
                        Equal(X, Times(Inverse(Times(A, B, C)), D)),
                        Equal(Y, Times(B, C))
                        )


class Example056(object):
    def __init__(self):

        # TAGS
        # CSE plus

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        D = Matrix("D", (n1, n1))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)
        C.set_property(properties.INPUT)
        D.set_property(properties.INPUT)

        X = Matrix("X", (n1, n1))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n1, n1))
        Y.set_property(properties.OUTPUT)

        Z = Matrix("Z", (n1, n1))
        Z.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(X, Plus(A, B, C,)),
                            Equal(Y, Plus(B, C, D)),
                            Equal(Z, Plus(D, A, C))
                            )


class Example057(object):
    def __init__(self):

        # TAGS
        # CSE plus

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        D = Matrix("D", (n1, n1))
        E = Matrix("E", (n1, n1))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)
        C.set_property(properties.INPUT)
        D.set_property(properties.INPUT)
        E.set_property(properties.INPUT)

        X = Matrix("X", (n1, n1))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n1, n1))
        Y.set_property(properties.OUTPUT)

        Z = Matrix("Z", (n1, n1))
        Z.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(X, Plus(A, B, C, Transpose(A), Transpose(B))),
                            # Equal(X, Plus(A, B, D, Transpose(A), Transpose(B), E, C)),
                            # Equal(Y, Plus(A, B, D, Transpose(A), Transpose(B))),
                            # Equal(Y, Plus(A, B, D, Transpose(A), Transpose(B), E)),
                            # Equal(Z, Plus(A, B)),
                            # Equal(Z, Plus(D, E)),
                            # Equal(X, Plus(C, Transpose(A), Transpose(B)))
                            )


class Example058(object):
    def __init__(self):

        # TAGS
        # CSE plus

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        D = Matrix("D", (n1, n1))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)
        C.set_property(properties.INPUT)
        D.set_property(properties.INPUT)
        D.set_property(properties.DIAGONAL)

        X = Matrix("X", (n1, n1))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n1, n1))
        Y.set_property(properties.OUTPUT)

        Z = Matrix("Z", (n1, n1))
        Z.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(X, Plus(A, B, D)),
                            Equal(Y, Plus(C, B, D)),
                            Equal(Z, Plus(A, B)),
                            )


class Example059(object):
    def __init__(self):

        # TAGS
        # SOP, POS

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        D = Matrix("D", (n1, n1))
        E = Matrix("E", (n1, n1))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)
        C.set_property(properties.INPUT)
        D.set_property(properties.INPUT)
        E.set_property(properties.INPUT)

        X = Matrix("X", (n1, n1))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n1, n1))
        Y.set_property(properties.OUTPUT)

        Z = Matrix("Z", (n1, n1))
        Z.set_property(properties.OUTPUT)

        # X = (A+D)B(C+E)
        self.eqns = Equations(
                            Equal(X, Times(Plus(A, D), B, Plus(C, E))),
                            # Equal(X, Plus(Times(A, B, C), Times(D, B, C), Times(D, B, E), Times(A, B, E))),
                            )


class Example060(object):
    def __init__(self):

        # TAGS
        # POS, SOP

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        D = Matrix("D", (n1, n1))
        E = Matrix("E", (n1, n1))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)
        C.set_property(properties.INPUT)
        D.set_property(properties.INPUT)
        E.set_property(properties.INPUT)

        X = Matrix("X", (n1, n1))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n1, n1))
        Y.set_property(properties.OUTPUT)

        Z = Matrix("Z", (n1, n1))
        Z.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(X, Plus(Times(A, B), Times(A, C))),
                            )


class Example061(object):
    def __init__(self):

        # TAGS
        # least squares, QR

        n = 10
        m = 5

        X = Matrix("X", (n, m))
        X.set_property(properties.INPUT)
        X.set_property(properties.FULL_RANK)
        X.set_property(properties.COLUMN_PANEL)

        y = Vector("y", (n, sONE))
        y.set_property(properties.INPUT)

        b = Vector("b", (m, sONE))
        b.set_property(properties.OUTPUT)

        # b = (X^T X)^-1 X^T y
        self.eqns = Equations(Equal(b, Times(Inverse(Times(Transpose(X), X)), Transpose(X), y)))


class Example062(object):
    def __init__(self):

        # TAGS
        # sum with scalars

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        alpha = Scalar("alpha")
        beta = Scalar("beta")

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)
        C.set_property(properties.INPUT)

        X = Matrix("X", (n1, n1))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n1, n1))
        Y.set_property(properties.OUTPUT)

        Z = Matrix("Z", (n1, n1))
        Z.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(X, Plus(Times(alpha, A), Times(beta, B), C)),
                            )


class Example063(object):
    def __init__(self):

        # TAGS
        # GWAS, generalized least squares

        n = 10 # 10
        m = 5 # 5

        S = Matrix("S", (n, n))
        S.set_property(properties.SPD)
        S.set_property(properties.INPUT)

        X = Matrix("X", (n, m))
        # X.set_property(properties.COLUMN_PANEL)
        X.set_property(properties.FULL_RANK)
        X.set_property(properties.INPUT)

        z = Vector("z", (m, sONE))
        z.set_property(properties.OUTPUT)

        y = Vector("y", (n, sONE))
        y.set_property(properties.INPUT)

        # clak.special_properties.add_expression(Times(Transpose(X), Inverse(S), X ), set([properties.SPD))

        self.eqns = Equations(Equal(z, Times(Inverse(Times(Transpose(X), Inverse(S), X ) ), Transpose(X), Inverse(S), y)))


class Example064(object):
    def __init__(self):

        # TAGS
        # transposed kernels

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.INPUT)

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)
        B.set_property(properties.INPUT)

        x = Vector("x", (n, sONE))
        x.set_property(properties.INPUT)

        y = Vector("y", (n, sONE))
        y.set_property(properties.INPUT)

        w = Vector("w", (sONE, n))
        w.set_property(properties.OUTPUT)

        self.eqns = Equations(
                        Equal(w, Times(Transpose(x), B)),
                        )


class Example065(object):
    def __init__(self):

        # TAGS
        # CSE times

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        D = Matrix("D", (n1, n1))
        E = Matrix("E", (n1, n1))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)
        C.set_property(properties.INPUT)
        D.set_property(properties.INPUT)
        E.set_property(properties.INPUT)

        X = Matrix("X", (n1, n1))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n1, n1))
        Y.set_property(properties.OUTPUT)

        Z = Matrix("Z", (n1, n1))
        Z.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(X, Times(A, B, C, D)),
                            Equal(Y, Times(B, C, D, E)),
                            )


class Example066(object):
    def __init__(self):

        # TAGS
        # CSE times

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        D = Matrix("D", (n1, n1))
        E = Matrix("E", (n1, n1))
        F = Matrix("F", (n1, n1))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)
        C.set_property(properties.INPUT)
        D.set_property(properties.INPUT)
        E.set_property(properties.INPUT)
        F.set_property(properties.INPUT)

        X = Matrix("X", (n1, n1))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n1, n1))
        Y.set_property(properties.OUTPUT)

        Z = Matrix("Z", (n1, n1))
        Z.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(X, Times(A, B, C, D)),
                            Equal(Y, Times(B, C, D, E)),
                            Equal(Z, Times(B, C, D, F)),
                            )


class Example067(object):
    def __init__(self):

        # TAGS
        # CSE times

        n1 = 1000
        n2 = 1000

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        D = Matrix("D", (n1, n1))
        E = Matrix("E", (n1, n1))
        F = Matrix("F", (n1, n1))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)
        C.set_property(properties.INPUT)
        D.set_property(properties.INPUT)
        E.set_property(properties.INPUT)
        F.set_property(properties.INPUT)

        X = Matrix("X", (n1, n1))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n1, n1))
        Y.set_property(properties.OUTPUT)

        Z = Matrix("Z", (n1, n1))
        Z.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(X, Times(A, B, C, D)),
                            Equal(Y, Times(B, C, D, E)),
                            Equal(Z, Times(B, C, D, E, F)),
                            )


class Example068(object):
    def __init__(self):

        # TAGS
        # matrix chain

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))

        A.set_property(properties.INPUT)
        A.set_property(properties.LOWER_TRIANGULAR)
        A.set_property(properties.SQUARE)
        B.set_property(properties.INPUT)
        B.set_property(properties.LOWER_TRIANGULAR)
        B.set_property(properties.SQUARE)
        C.set_property(properties.INPUT)

        X = Matrix("X", (n1, n1))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n1, n1))
        Y.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(Y, Times(C, InverseTranspose(A), InverseTranspose(B)))
                            )


class Example069(object):
    def __init__(self):

        # TAGS
        # CSE times

        n1 = 10
        n2 = 20

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))

        A.set_property(properties.INPUT)
        A.set_property(properties.LOWER_TRIANGULAR)
        A.set_property(properties.SQUARE)
        B.set_property(properties.INPUT)
        B.set_property(properties.LOWER_TRIANGULAR)
        B.set_property(properties.SQUARE)
        C.set_property(properties.INPUT)

        X = Matrix("X", (n1, n1))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n1, n1))
        Y.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(X, Times(Transpose(B), Transpose(A), C, Transpose(Inverse(A)), Transpose(Inverse(B)))),
                            )


class Example070(object):
    def __init__(self):

        # TAGS
        # complex inverse, random, explicit inverse

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.INPUT)

        D = Matrix("D", (n, n))
        D.set_property(properties.DIAGONAL)
        D.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        W.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(W, Times(Inverse(Plus(A, Inverse(D) ) ), D)))


class Example071(object):
    def __init__(self):

        # TAGS
        # complex inverse, random

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.INPUT)

        D = Matrix("D", (n, n))
        D.set_property(properties.DIAGONAL)
        D.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        W.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(W, Times(Inverse(Plus(A, D) ), D)))


class Example072(object):
    def __init__(self):

        # TAGS
        # complex inverse, termination/progress is a problem

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.INPUT)

        B = Matrix("B", (n, n))
        B.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        W.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(W, Plus(Inverse(A), B, A)))


class Example073(object):
    def __init__(self):

        # TAGS
        # explicit inverse, simple

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.INPUT)

        D = Matrix("D", (n, n))
        D.set_property(properties.DIAGONAL)
        D.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        W.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(W, Plus(Inverse(D), A)))


class Example074(object):
    def __init__(self):

        # TAGS
        # eigen trick

        n = 10

        Q = Matrix("Q", (n, n))
        Q.set_property(properties.ORTHOGONAL)
        Q.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        W.set_property(properties.INPUT)
        W.set_property(properties.DIAGONAL)

        I = IdentityMatrix(n, n)

        alpha = Scalar("alpha")

        two = ConstantScalar(2)

        X = Matrix("X", (n, n))
        X.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(X, Plus(Times(Transpose(Q), W, Q), Times(Plus(two, alpha), I))))


class Example075(object):
    def __init__(self):

        # TAGS
        # canonical form, simplify, SOP, POS, scalar in sum

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.INPUT)

        alpha = Scalar("alpha")

        two = ConstantScalar(2)

        X = Matrix("X", (n, n))
        X.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(X, Plus(Times(alpha, A), Times(two, A))))


class Example076(object):
    def __init__(self):

        # TAGS
        # eigen trick, SOP, POS

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.INPUT)

        B = Matrix("B", (n, n))
        B.set_property(properties.INPUT)

        Q = Matrix("Q", (n, n))
        Q.set_property(properties.ORTHOGONAL)
        Q.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        W.set_property(properties.INPUT)
        W.set_property(properties.DIAGONAL)

        I = IdentityMatrix(n, n)

        alpha = Scalar("alpha")

        two = ConstantScalar(2)

        X = Matrix("X", (n, n))
        X.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(X, Times(A, Plus(Times(Transpose(Q), W, Q), Times(alpha, I)), B)))


class Example077(object):
    def __init__(self):

        # TAGS
        # eigen trick, equivalent expression, SOP, POS

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.INPUT)

        B = Matrix("B", (n, n))
        B.set_property(properties.INPUT)

        Q = Matrix("Q", (n, n))
        Q.set_property(properties.ORTHOGONAL)
        Q.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        W.set_property(properties.INPUT)
        W.set_property(properties.DIAGONAL)

        I = IdentityMatrix(n, n)

        alpha = Scalar("alpha")

        two = ConstantScalar(2)

        X = Matrix("X", (n, n))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n, n))
        Y.set_property(properties.OUTPUT)

        self.eqns = Equations(
                        Equal(Y, Plus(A, B)),
                        Equal(X, Times(A, Plus(Times(Transpose(Q), W, Q), Times(alpha, I)), B))
                        )


class Example078(object):
    def __init__(self):

        # TAGS
        # symmetric product, CSE times

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        D = Matrix("D", (n1, n1))
        S = Matrix("S", (n1, n1))


        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)
        C.set_property(properties.INPUT)
        D.set_property(properties.INPUT)
        S.set_property(properties.INPUT)
        S.set_property(properties.SPD)

        X = Matrix("X", (n1, n1))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n1, n1))
        Y.set_property(properties.OUTPUT)

        Z = Matrix("Z", (n1, n1))
        Z.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(X, Times(A, B, S, Transpose(B), Transpose(A))),
                            )


class Example079(object):
    def __init__(self):

        # TAGS
        # CSE times

        n1 = 10
        n2 = 20

        A = Matrix("A", (n2, n1))
        B = Matrix("B", (n1, n1))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)

        X = Matrix("X", (n2, n1))
        X.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(X, Times(Inverse(Times(A, B, Transpose(B), Transpose(A))), A, B)),
                            )


class Example080(object):
    def __init__(self):

        # TAGS
        # CSE times

        n1 = 10
        n2 = 20

        A = Matrix("A", (n2, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n2, n2))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)

        X = Matrix("X", (n2, n1))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n2, n2))
        Y.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(X, Times(C, A, B)),
                            Equal(Y, Times(A, B, Transpose(B), Transpose(A)))
                            )


class Example081(object):
    def __init__(self):

        # TAGS
        # CSE times

        n = 10

        A = Matrix("A", (n, n))

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)
        B.set_property(properties.INPUT)

        C = Matrix("C", (n, n))
        C.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        W.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(W, Times(InverseTranspose(B), InverseTranspose(C), A, Transpose(C), Transpose(B))))


class Example082(object):
    def __init__(self):

        # TAGS
        # CSE times

        n = 10

        A = Matrix("A", (n, n))

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)
        B.set_property(properties.INPUT)

        C = Matrix("C", (n, n))
        C.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        W.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(W, Times(Transpose(B), InverseTranspose(C), A, Transpose(C), InverseTranspose(B))))


class Example083(object):
    def __init__(self):

        # TAGS
        # CSE times

        n = 10

        A = Matrix("A", (n, n))
        D = Matrix("D", (n, n))

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)
        B.set_property(properties.INPUT)

        C = Matrix("C", (n, n))
        C.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        W.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(W, Times(InverseTranspose(B), InverseTranspose(C), A, Transpose(C), Transpose(B),  D, Inverse(C), Inverse(B))))


class Example084(object):
    def __init__(self):

        # TAGS
        # CSE times, termination/progress is a problem

        n = 10

        A = Matrix("A", (n, n))

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)
        B.set_property(properties.INPUT)

        C = Matrix("C", (n, n))
        C.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        W.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(W, Plus(Times(B, C), Times(Inverse(B), A))))


class Example085(object):
    def __init__(self):

        # TAGS
        # CSE times, symmetric product

        n1 = 10
        n2 = 10

        A = Matrix("A", (n2, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n2, n2))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)

        X = Matrix("X", (n2, n1))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n1, n1))
        Y.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(Y, Times(Transpose(A), Transpose(B), B, A))
                            )

class Example086(object):
    def __init__(self):

        # TAGS
        # CSE times, symmetric product

        n1 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)

        A.set_property(properties.LOWER_TRIANGULAR)
        B.set_property(properties.LOWER_TRIANGULAR)
        C.set_property(properties.LOWER_TRIANGULAR)

        X = Matrix("X", (n1, n1))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n1, n1))
        Y.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(X, Plus(Times(A, B, Transpose(B), Transpose(A)), C))
                            )

class Example087(object):
    def __init__(self):

        # TAGS
        # matrix chain, indices, example for the matrix chain paper

        i = Index("i", 2)
        j = Index("j", 2)
        n1 = 10000
        n2 = 11
        n3 = 10
        n4 = 10
        n5 = 10

        A = Matrix("A", (n1, n2), set())
        B = Matrix("B", (n2, n3), set())
        C = Matrix("C", (n3, n4), set([i]))
        D = Matrix("D", (n4, n5), set([j]))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)
        C.set_property(properties.INPUT)
        D.set_property(properties.INPUT)

        X = Matrix("X", (n1, n5), set([i, j]))

        X.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(X, Times(A, B, C, D)) 
                            )

class Example088(object):
    def __init__(self):

        # TAGS
        # trick, tricks

        n1 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)
        C.set_property(properties.INPUT)

        X = Matrix("X", (n1, n1))

        X.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(X, Plus(Times(Transpose(A), B), Times(Transpose(B), A), Times(Transpose(A), A), C)) 
                            )

class Example089(object):
    def __init__(self):

        # TAGS
        # eigen trick

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.INPUT)

        B = Matrix("B", (n, n))
        B.set_property(properties.INPUT)

        Q = Matrix("Q", (n, n))
        Q.set_property(properties.ORTHOGONAL)
        Q.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        W.set_property(properties.INPUT)
        W.set_property(properties.DIAGONAL)

        I = IdentityMatrix(n, n)

        alpha = Scalar("alpha")

        two = ConstantScalar(2)

        X = Matrix("X", (n, n))
        X.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(X, Plus(Times(Transpose(Q), W, Q), Times(alpha, I), B)))


class Example090(object):
    def __init__(self):

        # TAGS
        # blocks, blocked, partitioning

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.INPUT)

        B = Matrix("B", (n, n))
        B.partitioning = (set([5]), set([5]))
        B.set_property(properties.INPUT)

        D = Matrix("D", (n, n))
        D.set_property(properties.INPUT)
        # D.set_property(properties.DIAGONAL)

        alpha = Scalar("alpha")

        I = IdentityMatrix(n, n)

        X = Matrix("X", (n, n))
        X.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(X, Times(alpha, B, D)))


class Example091(object):
    def __init__(self):

        # TAGS
        # blocks, blocked, partitioning, indices

        n = 10

        i = Index("i", 50)

        A = Matrix("A", (n, n))
        A.set_property(properties.INPUT)

        B = Matrix("B", (n, n), set([i]))
        B.set_property(properties.INPUT)

        C = Matrix("C", (2*n, n))
        C.set_property(properties.INPUT)
        # D.set_property(properties.DIAGONAL)

        alpha = Scalar("alpha")

        I = IdentityMatrix(n, n)

        X = Matrix("X", (n, n), set([i]))
        X.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(X, Times(BlockedExpression([[A, B]]), C)))


class Example092(object):
    def __init__(self):

        # TAGS
        # explicit inverse, matrix chain

        n = 100
        m = 10

        V = Matrix("V", (n, m))
        V.set_property(properties.INPUT)
        V.set_property(properties.FULL_RANK)

        X = Matrix("X", (m, m))
        X.set_property(properties.OUTPUT)


        self.eqns = Equations(Equal(X, Inverse(Times(Transpose(V), V))))
        # self.eqns = Equations(Equal(X, Times(Transpose(V), V)))


class Example093(object):
    def __init__(self):

        # TAGS
        # CSE times, symmetric

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.ORTHOGONAL)

        B = Matrix("B", (n, n))
        B.set_property(properties.DIAGONAL)

        C = Matrix("C", (n, n))
        D = Matrix("D", (n, n))

        X = Matrix("X", (n, n))
        X.set_property(properties.OUTPUT)


        self.eqns = Equations(Equal(X, Times(A, Inverse(B), Transpose(A), C)),
                              Equal(X, Times(D, A, Inverse(B), Transpose(A)))
                              )


class Example094(object):
    def __init__(self):

        # TAGS
        # indexed matrix chain, indices

        i = Index("i", 10)
        j = Index("j", 20)
        k = Index("k", 100)
        n1 = 600
        n2 = 600
        n3 = 100
        n4 = 500
        n5 = 200
        n6 = 10
        n7 = 10

        A = Matrix("A", (n1, n2), set([i]))
        B = Matrix("B", (n2, n3), set([i, j]))
        C = Matrix("C", (n3, n4), set())
        D = Matrix("D", (n4, n5), set())
        # E = Matrix("E", (n5, n6), set([i, j, k]))
        # F = Matrix("F", (n6, n7), set())
        # # G = Matrix("G", (10, 10), set([]))
        # # H = Matrix("H", (10, 10), set([]))

        A.set_property(properties.INPUT)
        A.set_property(properties.SPD)
        B.set_property(properties.INPUT)
        C.set_property(properties.INPUT)
        D.set_property(properties.INPUT)
        # E.set_property(properties.INPUT)
        # F.set_property(properties.INPUT)
        # G.set_property(properties.INPUT)
        # H.set_property(properties.INPUT)

        X = Matrix("X", (n1, n5), set([i, j]))

        X.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            # Equal(X, Times(A, B, C, D, E, F, G, H))
                            Equal(X, Times(Inverse(A), B, C, D)) 
                            )

class Example095(object):
    def __init__(self):

        # TAGS
        # indexed matrix chain, indices

        i = Index("i", 2)
        j = Index("j", 20)
        n1 = 10
        n2 = 10
        n3 = 10

        A = Matrix("A", (n1, n2), set([i]))
        B = Matrix("B", (n2, n3), set())
        c = Vector("c", (n3, 1), set([j]))


        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)
        c.set_property(properties.INPUT)

        X = Matrix("X", (n1, 1), set([i, j]))

        X.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(X, Times(A, B, c)) 
                            )

class Example096(object):
    def __init__(self):

        # TAGS
        # simple, matrix chain

        n1 = 1500
        n2 = 1000
        n3 = 500

        A = Matrix("A", (n1, n2))
        A.set_property(properties.INPUT)

        B = Matrix("B", (n2, n2))
        B.set_property(properties.LOWER_TRIANGULAR)
        B.set_property(properties.INPUT)

        C = Matrix("C", (n2, n2))
        C.set_property(properties.UPPER_TRIANGULAR)
        C.set_property(properties.INPUT)

        D = Matrix("D", (n2, n3))
        D.set_property(properties.DIAGONAL)
        D.set_property(properties.INPUT)

        W = Matrix("W", (n1, n3))
        W.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(W, Times(A, Inverse(B), InverseTranspose(C), D)))


class Example097(object):
    def __init__(self):

        # TAGS
        # simple, matrix chain

        n1 = 600
        n2 = 200
        n3 = 1500
        n4 = 800
        n5 = 1000

        A = Matrix("A", (n1, n1))
        A.set_property(properties.LOWER_TRIANGULAR)
        A.set_property(properties.INPUT)

        B = Matrix("B", (n1, n2))
        B.set_property(properties.INPUT)

        C = Matrix("C", (n2, n3))
        C.set_property(properties.INPUT)

        D = Matrix("D", (n3, n3))
        D.set_property(properties.UPPER_TRIANGULAR)
        D.set_property(properties.INPUT)

        E = Matrix("E", (n3, n4))
        E.set_property(properties.UPPER_TRIANGULAR)
        E.set_property(properties.INPUT)

        F = Matrix("F", (n4, n5))
        F.set_property(properties.INPUT)

        W = Matrix("W", (n1, n5))
        W.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(W, Times(Inverse(A), B, C, InverseTranspose(D), E, F)))


class Example098(object):
    def __init__(self):

        # TAGS
        # simple, scalar, constant

        alpha = Scalar("alpha")
        two = ConstantScalar(2)

        beta = Scalar("beta")

        self.eqns = Equations(Equal(beta, Times(alpha, two)))


class Example099(object):
    def __init__(self):

        # TAGS
        # simple, scalar, constant

        alpha = Scalar("alpha")
        two = ConstantScalar(2)

        beta = Scalar("beta")

        self.eqns = Equations(Equal(beta, Plus(alpha, two)))


class Example100(object):
    def __init__(self):

        # TAGS
        # simple, scalar, constant

        n = 10

        A = Matrix("A", (n, n))
        two = ConstantScalar(2)

        X = Matrix("X", (n, n))


        beta = Scalar("beta")

        self.eqns = Equations(Equal(X, Times(two, A)))


class Example101(object):
    def __init__(self):

        # TAGS
        # 

        n1 = 100
        n2 = 1000
        n3 = 300

        A = Matrix("A", (n2, n1))
        B = Matrix("B", (n3, n1))


        X = Matrix("X", (n1, n1))


        self.eqns = Equations(Equal(X, Plus(Times(Transpose(A), A), Times(Transpose(B), B))))
        # self.eqns = Equations(Equal(X, Plus(Times(Transpose(A), A), Times(Transpose(B), B), IdentityMatrix(n1, n1))))


class Example102(object):
    def __init__(self):

        # TAGS
        # 

        n1 = 100
        n2 = 1000
        n3 = 300

        A = Matrix("A", (n2, n1))
        B = Matrix("B", (n3, n1))


        X = Matrix("X", (n1, n1))


        self.eqns = Equations(Equal(X, Times(Transpose(A), A, Transpose(B), B)))


class Example103(object):
    def __init__(self):

        # TAGS
        # 

        n1 = 50
        n2 = 1000
        n3 = 250

        A = Matrix("A", (n2, n1))
        A.set_property(properties.FULL_RANK)
        B = Matrix("B", (n3, n1))
        B.set_property(properties.FULL_RANK)
        y = Vector("y", (n1, 1))

        x = Vector("x", (n1, 1))

        self.eqns = Equations(Equal(x, Times(Inverse(Times(Transpose(A), A, Transpose(B), B)), y)))


class Example104(object):
    def __init__(self):

        # TAGS
        # CSE general

        n1 = 50
        n2 = 60
        n3 = 70

        A = Matrix("A", (n1, n2))
        B = Matrix("B", (n2, n3))
        X = Matrix("X", (n3, n1))

        Y = Matrix("Y", (n3, n1))

        self.eqns = Equations(Equal(X, Inverse(Times(A, B))),
                              Equal(Y, Inverse(Times(A, B))))



class Example105(object):
    def __init__(self):

        # TAGS
        # simple, plus, scalars

        n1 = 50

        alpha = Scalar("alpha")
        beta = Scalar("beta")
        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        X = Matrix("X", (n1, n1))


        self.eqns = Equations(Equal(X, Plus(Times(alpha, A), Times(beta, B))))


class Example106(object):
    def __init__(self):

        # TAGS
        # simple, times

        n1 = 50

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        X = Matrix("X", (n1, n1))

        # self.eqns = Equations(Equal(X, Plus(A, B, C)))
        self.eqns = Equations(Equal(X, Plus(Transpose(A), Transpose(B))))


class Example107(object):
    def __init__(self):

        # TAGS
        # simple, plus, identity

        n1 = 50

        A = Matrix("A", (n1, n1))
        K = Matrix("K", (n1, n1))
        X = Matrix("X", (n1, n1))


        self.eqns = Equations(Equal(X, Plus(A, IdentityMatrix(n1, n1))))


class Example108(object):
    def __init__(self):

        # TAGS
        # simple, plus, scalars

        n1 = 50

        alpha = Scalar("alpha")
        beta = Scalar("beta")
        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        X = Matrix("X", (n1, n1))


        self.eqns = Equations(Equal(X, Plus(A, Times(Plus(alpha, beta), B))))


class Example109(object):
    def __init__(self):

        # TAGS
        # simple, plus, scalars

        n1 = 50

        alpha = Scalar("alpha")
        beta = Scalar("beta")
        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        X = Matrix("X", (n1, n1))


        self.eqns = Equations(Equal(X, Plus(A, Times(alpha, beta, B))))


class Example110(object):
    def __init__(self):

        # TAGS
        # eigen trick

        n = 10

        Q = Matrix("Q", (n, n))
        Q.set_property(properties.ORTHOGONAL)
        Q.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        W.set_property(properties.INPUT)
        W.set_property(properties.DIAGONAL)

        I = IdentityMatrix(n, n)

        alpha = Scalar("alpha")

        two = ConstantScalar(2)

        X = Matrix("X", (n, n))
        X.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(X, Plus(Times(Transpose(Q), W, Q), Times(alpha, I))))


class Example111(object):
    def __init__(self):

        # TAGS
        # 

        n1 = 50
        n2 = 1000
        n3 = 250

        A = Matrix("A", (n2, n1))
        A.set_property(properties.FULL_RANK)
        B = Matrix("B", (n3, n1))
        B.set_property(properties.FULL_RANK)
        y = Vector("y", (n1, 1))

        x = Vector("x", (n1, 1))

        # x = inv()
        self.eqns = Equations(Equal(x, Times(Inverse(Plus(Times(Transpose(A), A), Times(Transpose(B), B))), y)))


class Example112(object):
    def __init__(self):

        # TAGS
        # eigen trick

        n1 = 50

        alpha = Scalar("alpha")

        I = IdentityMatrix(n1, n1)

        S = Matrix("S", (n1, n1))
        S.set_property(properties.SPD)

        y = Vector("y", (n1, 1))

        x = Vector("x", (n1, 1))


        self.eqns = Equations(Equal(x, Times(Inverse(Plus(S, Times(alpha, I))), y)))


class Example113(object):
    def __init__(self):

        # TAGS
        # POS, SOP

        n1 = 10
        n2 = 10 

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        D = Matrix("D", (n1, n1))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)
        C.set_property(properties.INPUT)
        D.set_property(properties.INPUT)

        X = Matrix("X", (n1, n1))
        X.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(X, Times(Plus(A, B), C))
                            )


class Example114(object):
    def __init__(self):

        # TAGS
        # simple, product

        n1 = 10
        n2 = 10 

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))

        X = Matrix("X", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Times(A, B))
                            )


class Example115(object):
    def __init__(self):

        # TAGS
        # simple, matrix chain

        n1 = 600
        n2 = 200
        n3 = 1500
        n4 = 800
        n5 = 1000

        A = Matrix("A", (n1, n1))
        A.set_property(properties.LOWER_TRIANGULAR)
        A.set_property(properties.INPUT)

        B = Matrix("B", (n1, n2))
        B.set_property(properties.INPUT)

        C = Matrix("C", (n2, n3))
        C.set_property(properties.INPUT)

        D = Matrix("D", (n3, n3))
        D.set_property(properties.UPPER_TRIANGULAR)
        D.set_property(properties.INPUT)

        E = Matrix("E", (n3, n4))
        E.set_property(properties.UPPER_TRIANGULAR)
        E.set_property(properties.INPUT)

        F = Matrix("F", (n4, n5))
        F.set_property(properties.INPUT)

        W = Matrix("W", (n1, n5))
        W.set_property(properties.OUTPUT)

        self.eqns = Equations(Equal(W, Times(A, B, C, D, E, F)))


class Example116(object):
    def __init__(self):

        # TAGS
        # CSE plus

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        D = Matrix("D", (n1, n1))
        E = Matrix("E", (n1, n1))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)
        C.set_property(properties.INPUT)
        D.set_property(properties.INPUT)
        E.set_property(properties.INPUT)

        X = Matrix("X", (n1, n1))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n1, n1))
        Y.set_property(properties.OUTPUT)

        Z = Matrix("Z", (n1, n1))
        Z.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(X, Times(A, Transpose(B), C, B, Transpose(A))),
                            )


class Example117(object):
    def __init__(self):

        # TAGS
        # matrix chain, linear system

        n = 10

        B = Matrix("B", (n, n))
        B.set_property(properties.SPD)

        x = Vector("x", (n, sONE))

        w = Vector("w", (n, sONE))

        self.eqns = Equations(Equal(w, Times(Inverse(B), x)))


class Example118(object):
    def __init__(self):

        # TAGS
        # simple, sum, inverse

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.NON_SINGULAR)

        X = Matrix("X", (n, n))

        self.eqns = Equations(Equal(X, Plus(A, Inverse(A))))


class Example119(object):
    def __init__(self):

        # TAGS
        # simple, sum, inverse

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.SPD)

        X = Matrix("X", (n, n))

        self.eqns = Equations(Equal(X, Plus(A, Inverse(A))))


class Example120(object):
    def __init__(self):

        # TAGS
        # trick, tricks

        n1 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)
        C.set_property(properties.INPUT)

        X = Matrix("X", (n1, n1))

        X.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(X, Plus(Times(Transpose(A), B), Times(Transpose(B), A), Times(Transpose(A), A))) 
                            )


class Example121(object):
    def __init__(self):

        # TAGS
        # conjugate prior

        n1 = 10

        S = Matrix("S", (n1, n1))
        S.set_property(properties.SPD)

        P = Matrix("P", (n1, n1))
        P.set_property(properties.SPD)

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

class Example124(object):
    def __init__(self):

        # TAGS
        # CSE times

        n1 = 1000
        n2 = 1000

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        D = Matrix("D", (n1, n1))

        A.set_property(properties.INPUT)
        B.set_property(properties.INPUT)
        B.set_property(properties.NON_SINGULAR)
        # B.set_property(properties.SPD)
        B.set_property(properties.LOWER_TRIANGULAR)
        C.set_property(properties.INPUT)
        D.set_property(properties.INPUT)

        X = Matrix("X", (n1, n1))
        X.set_property(properties.OUTPUT)

        Y = Matrix("Y", (n1, n1))
        Y.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(X, Times(A, InverseTranspose(B), C)),
                            Equal(Y, Times(D, Inverse(B), Transpose(A))),
                            )


class Example125(object):
    def __init__(self):

        # TAGS
        # nino2016, application

        nsd = 1000
        msd = 1000
        N = 1000
        m = 1000
        
        minusone = ConstantScalar(-1)

        B = Matrix("B", (N, N))
        H = Matrix("H", (msd, N))
        R = Matrix("R", (msd, m))
        Y = Matrix("Y", (msd, N))
        Z = Matrix("Z", (nsd, N)) # originally X^b
        X = Matrix("X", (nsd, N))

        B.set_property(properties.INPUT)
        H.set_property(properties.INPUT)
        R.set_property(properties.INPUT)
        Y.set_property(properties.INPUT)
        Z.set_property(properties.INPUT)
        X.set_property(properties.OUTPUT)

        self.eqns = Equations(
                            Equal(
                                X,
                                Plus(
                                    Z,
                                    Times(
                                        Inverse(Plus(Inverse(B), Times(Transpose(H), Inverse(R), H))),
                                        Plus(Y, Times(minusone, H, Z))
                                        )
                                    )
                                )
                            )

        # self.eqns = Equations(
        #                     Equal(
        #                         X,
        #                         Times(
        #                                 Inverse(Plus(Inverse(B), Times(Transpose(H), Inverse(R), H))),
        #                                 Plus(Y, Times(H, Z))
        #                                 )
        #                             )
        #                     )
        # print(self.eqns)


class Example126(object):
    def __init__(self):

        # TAGS
        # straszak2015, application

        # n > m
        m = 1000
        n = 1500

        A = Matrix("A", (m, n))
        A.set_property(properties.FULL_RANK)
        A.set_property(properties.INPUT)

        W = Matrix("W", (n, n))
        # C = Matrix("C", (n1, n1))

        # W is positive
        W.set_property(properties.FULL_RANK)
        W.set_property(properties.DIAGONAL)
        W.set_property(properties.SPD)
        W.set_property(properties.SYMMETRIC)
        W.set_property(properties.NON_SINGULAR)
        W.set_property(properties.INPUT)

        b = Vector("b", (m, sONE))
        b.set_property(properties.INPUT)

        c = Vector("c", (n, sONE))
        c.set_property(properties.INPUT)

        x = Vector("x", (n, sONE))
        x.set_property(properties.OUTPUT)

        minusone = ConstantScalar(-1.0)

        # derivation.special_properties.add_expression(Times(A, W, Transpose(A)), {properties.SPD})

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