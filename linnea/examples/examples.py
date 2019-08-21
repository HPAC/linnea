

from ..algebra.expression import Symbol, Scalar, Vector, Matrix, ConstantScalar, \
                                 Equal, Plus, Times, Transpose, Inverse, \
                                 InverseTranspose, InverseConjugate, \
                                 InverseConjugateTranspose, \
                                 ConjugateTranspose, Index, IdentityMatrix

from ..algebra.properties import Property as properties

from ..algebra.equations import Equations

from .. import derivation


#######################
### Diego's examples

class Example001():
    def __init__(self):

        # TAGS
        # matrix chain

        n = 10

        x = Vector("x", (n, 1))

        y = Vector("y", (n, 1))

        z = Vector("z", (n, 1))

        alpha = Scalar("alpha")

        # alpha = x^T z x^T y
        self.eqns = Equations(Equal(alpha, Times(Transpose(x), z, Transpose(x), y)))

        

###

class Example002():
    def __init__(self):

        # TAGS
        # SOP, POS, CSE

        n = 10

        X = Matrix("X", (n, n))

        Y = Matrix("Y", (n, n))

        v = Vector("v", (n, 1))

        w = Vector("w", (n, 1))

        # w = X Y^T v + Y X^T v
        self.eqns = Equations(Equal(w, Plus(Times(X, Transpose(Y), v), Times(Y, Transpose(X), v))))

###

class Example003():
    def __init__(self):

        # TAGS
        # matrix chain

        n = 10

        L = Matrix("L", (n, n))
        L.set_property(properties.LOWER_TRIANGULAR)

        v = Vector("v", (n, 1))

        u = Vector("u", (n, 1))

        alpha = Scalar("alpha")

        # alpha = u^T L^-1 L^-1 v
        self.eqns = Equations(Equal(alpha, Times(Transpose(u), Inverse(L), Inverse(L), v)))

        expr = Times(Transpose(u), Inverse(L), Inverse(L), v)

###

class Example004():
    def __init__(self):

        # TAGS
        # matrix chain, CSE

        n = 10

        x = Vector("x", (n, 1))

        y = Vector("y", (n, 1))

        alpha = Scalar("alpha")

        # alpha = x^T y x^T y
        self.eqns = Equations(Equal(alpha, Times(Transpose(x), y, Transpose(x), y)))

###

class Example005():
    def __init__(self):

        # TAGS
        # solve equation

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.SPD)
        A.set_property(properties.SQUARE)
        A.set_property(properties.NON_SINGULAR)

        B = Matrix("B", (n, n))

        X = Matrix("X4", (n, n))

        # A X = B
        self.eqns = Equations(Equal(Times(A, X), B))

###

class Example006():
    def __init__(self):

        # TAGS
        # explicit inverse, matrix chain

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.SYMMETRIC)
        A.set_property(properties.SQUARE)
        A.set_property(properties.NON_SINGULAR)

        B = Matrix("B", (n, n))
        B.set_property(properties.SYMMETRIC)
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)

        # B = A^-1
        self.eqns = Equations(Equal(B, Inverse(A)))

###

class Example007():
    def __init__(self):

        # TAGS
        # least squares, QR

        n = 1500
        m = 1000

        X = Matrix("X", (n, m))
        X.set_property(properties.FULL_RANK)

        y = Vector("y", (n, 1))

        b = Vector("b", (m, 1))

        # b = (X^T X)^-1 X^T y
        self.eqns = Equations(Equal(b, Times(Inverse(Times(Transpose(X), X)), Transpose(X), y)))

###

class Example008():
    def __init__(self):

        # TAGS
        # simple

        n = 10

        A = Matrix("A", (n, n))

        B = Matrix("B", (n, n))

        C = Matrix("C", (n, n))

        D = Matrix("D", (n, n))

        E = Matrix("E", (n, n))

        X = Matrix("X", (n, n))

        # X = A B + C D + E
        self.eqns = Equations(Equal(X, Plus(Times(A, B), Times(C, D), E)))

###

class Example009():
    def __init__(self):

        # TAGS
        # CSE, matrix chain

        n = 10

        L = Matrix("L", (n, n))
        L.set_property(properties.LOWER_TRIANGULAR)

        A = Matrix("A", (n, n))

        X = Matrix("X", (n, n))

        # X = A^T L^-1 L^-T A
        self.eqns = Equations(Equal(X, Times(Transpose(A), Inverse(L), Transpose(Inverse(L)), A)))

### GWAS in Demo-Rev missing and derivative

###############################
### my examples


class Example010():
    def __init__(self):

        # TAGS
        # complex inverse

        n = 10
        m = 5

        S = Matrix("S", (n, n))
        S.set_property(properties.SPD)

        CP = Matrix("CP", (n, m))
        CP.set_property(properties.COLUMN_PANEL)
        CP.set_property(properties.FULL_RANK)

        Y = Matrix("Y", (m, m))

        self.eqns = Equations(Equal(Y, Inverse(Times(Transpose(CP), S, CP ) )))


class Example011():
    def __init__(self):

        # TAGS
        # symmetric product

        n = 10
        m = 5

        S = Matrix("S", (n, n))
        S.set_property(properties.SPD)

        CP = Matrix("CP", (n, m))
        CP.set_property(properties.COLUMN_PANEL)
        CP.set_property(properties.FULL_RANK)

        Y = Matrix("Y", (m, m))

        self.eqns = Equations(Equal(Y, Times(Transpose(CP), S, CP )))


class Example012():
    def __init__(self):
        
        # TAGS
        # simple

        n = 10

        A = Matrix("A", (n, n))

        W = Matrix("W", (n, n))

        self.eqns = Equations(Equal(W, Times(Transpose(A), A)))


class Example013():
    def __init__(self):

        # TAGS
        # simple, matrix chain

        n = 10

        A = Matrix("A", (n, n))

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)

        C = Matrix("C", (n, n))

        D = Matrix("D", (n, n))
        D.set_property(properties.DIAGONAL)

        W = Matrix("W", (n, n))

        self.eqns = Equations(Equal(W, Times(A, B, C, D)))

class Example014():
    def __init__(self):

        # TAGS
        # matrix chain, CSE, symmetric product

        n = 10

        A = Matrix("A", (n, n))

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)

        W = Matrix("W", (n, n))

        self.eqns = Equations(Equal(W, Times(Transpose(B), Transpose(A), A, B)))


class Example015():
    def __init__(self):

        # TAGS
        # simple

        n = 10

        A = Matrix("A", (n, n))

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)

        C = Matrix("C", (n, n))

        D = Matrix("D", (n, n))
        D.set_property(properties.DIAGONAL)

        W = Matrix("W", (n, n))

        self.eqns = Equations(Equal(W, Plus(A, B, C, D)))


class Example016():
    def __init__(self):

        # TAGS
        # complex inverse, random

        n = 10

        A = Matrix("A", (n, n))

        D = Matrix("D", (n, n))
        D.set_property(properties.DIAGONAL)

        S = Matrix("S", (n, n))
        S.set_property(properties.SPD)

        W = Matrix("W", (n, n))

        self.eqns = Equations(Equal(W, Times(A, Inverse(S), Inverse(Plus(A, Inverse(S) ) ), D)))


class Example017():
    def __init__(self):

        # TAGS
        # complex inverse, random

        n = 10

        A = Matrix("A", (n, n))

        D = Matrix("D", (n, n))
        D.set_property(properties.DIAGONAL)

        S = Matrix("S", (n, n))
        S.set_property(properties.SPD)

        W = Matrix("W", (n, n))

        self.eqns = Equations(Equal(W, Times(A, Inverse(S), Inverse(Plus(A, Inverse(D) ) ), D)))


class Example018():
    def __init__(self):

        # TAGS
        # random, explicit inverse, sum

        n = 10

        A = Matrix("A", (n, n))

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)

        C = Matrix("C", (n, n))

        D = Matrix("D", (n, n))
        D.set_property(properties.DIAGONAL)

        W = Matrix("W", (n, n))

        self.eqns = Equations(Equal(W, Plus(A, Inverse(B), Inverse(C), D)))


class Example019():
    def __init__(self):

        # TAGS
        # random, explicit inverse, sum

        n = 10

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)

        C = Matrix("C", (n, n))

        D = Matrix("D", (n, n))
        D.set_property(properties.DIAGONAL)

        W = Matrix("W", (n, n))

        self.eqns = Equations(Equal(W, Plus(Inverse(B), Inverse(C), D)))


class Example020():
    def __init__(self):

        # TAGS
        # explicit inverse

        n = 10

        S = Matrix("S", (n, n))
        S.set_property(properties.SPD)

        W = Matrix("W", (n, n))

        self.eqns = Equations(Equal(W, Inverse(S)))


class Example021():
    def __init__(self):

        # TAGS
        # matrix chain, linear system

        n = 10

        A = Matrix("A", (n, n))

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)
        # B.set_property(properties.SPD)

        x = Vector("x", (n, 1))

        w = Vector("w", (n, 1))

        self.eqns = Equations(Equal(w, Times(A, Inverse(B), x)))


class Example022():
    def __init__(self):

        # TAGS
        # linear system, matrix chain

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.SQUARE)
        A.set_property(properties.NON_SINGULAR)

        x = Vector("x", (n, 1))

        w = Vector("w", (n, 1))

        self.eqns = Equations(Equal(w, Times(Inverse(A), x)))


class Example023():
    def __init__(self):

        # TAGS
        # triangular system

        n = 1000

        L = Matrix("L", (n, n))
        L.set_property(properties.LOWER_TRIANGULAR)
        L.set_property(properties.NON_SINGULAR)

        x = Vector("x", (n, 1))

        w = Vector("w", (n, 1))

        self.eqns = Equations(Equal(w, Times(Inverse(L), x)))


class Example024():
    def __init__(self):

        # TAGS
        # SPD system

        n = 10

        S = Matrix("S", (n, n))
        S.set_property(properties.SPD)

        x = Vector("x", (n, 1))

        w = Vector("w", (n, 1))

        self.eqns = Equations(Equal(w, Times(Inverse(S), x)))


class Example025():
    def __init__(self):

        # TAGS
        # CSE, SOP, POS, random

        n = 10

        A = Matrix("A", (n, n))

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)

        S = Matrix("S", (n, n))
        S.set_property(properties.SPD)

        W = Matrix("W", (n, n))

        self.eqns = Equations(Equal(W, Plus(Times(A, B, S), Times(Transpose(B), Transpose(A), S) )))


class Example026():
    def __init__(self):

        # TAGS
        # random, sum, simple

        n = 10

        A = Matrix("A", (n, n))

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)

        C = Matrix("C", (n, n))

        D = Matrix("D", (n, n))
        D.set_property(properties.DIAGONAL)

        W = Matrix("W", (n, n))

        self.eqns = Equations(Equal(W, Plus(A, Transpose(B), C, D)))


class Example027():
    def __init__(self):

        # TAGS
        # transposed kernel

        n = 10

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)

        C = Matrix("C", (n, n))

        W = Matrix("W", (n, n))

        self.eqns = Equations(Equal(W, Times(Transpose(B), Transpose(C))))


class Example028():
    def __init__(self):

        # TAGS
        # transposed kernel, CSE

        n = 10

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)

        C = Matrix("C", (n, n))

        W = Matrix("W", (n, n))

        self.eqns = Equations(Equal(W, Times(Transpose(B), Transpose(C), C, B)))


class Example029():
    def __init__(self):

        # TAGS
        # transposed kernel

        n = 10

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)

        x = Vector("x", (n, 1))

        y = Vector("y", (n, 1))

        W = Matrix("W", (n, n))

        self.eqns = Equations(Equal(W, Times(y, Transpose(x), B)))


class Example030():
    def __init__(self):

        # TAGS
        # transposed kernel

        n = 10

        A = Matrix("A", (n, n))

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)

        x = Vector("x", (n, 1))

        y = Vector("y", (n, 1))

        W = Matrix("W", (n, n))

        U = Matrix("U", (n, n))

        self.eqns = Equations(
                        Equal(W, Times(y, Transpose(x), B)),
                        Equal(U, Times(Transpose(A), A))
                        )

class Example031():
    def __init__(self):

        # TAGS
        # blocked algorithm

        n = 10
        m = 10
        k = 10

        L00 = Matrix("L00", (n, n))
        L00.set_property(properties.LOWER_TRIANGULAR)

        L11 = Matrix("L11", (m, m))
        L11.set_property(properties.LOWER_TRIANGULAR)

        L22 = Matrix("L22", (k, k))
        L22.set_property(properties.LOWER_TRIANGULAR)

        L10 = Matrix("L10", (m, n))

        L21 = Matrix("L21", (k, m))

        X10_in = Matrix("X10_in", (m, n))

        X10 = Matrix("X10", (m, n))

        X11 = Matrix("X11", (m, m))

        X20_in = Matrix("X20_in", (k, n))

        X20 = Matrix("X20", (k, n))

        X21 = Matrix("X21", (k, m))

        self.eqns = Equations(
                        Equal(X10, Times(X10_in, Inverse(L00))),
                        Equal(X20, Plus(X20_in, Times(Inverse(L22), L21, Inverse(L11), L10))),
                        Equal(X11, Inverse(L11)),
                        # Equal(X21, Minus([Times(Inverse(L22), L21))), # this is the correct equation, but minus is a problem
                        Equal(X21, Times(Inverse(L22), L21)),
                        )
        

class Example032():
    def __init__(self):

        # TAGS
        # CSE times

        n = 10

        A = Matrix("A", (n, n))

        B = Matrix("B", (n, n))

        C = Matrix("C", (n, n))

        D = Matrix("D", (n, n))

        X = Matrix("X", (n, n))

        Y = Matrix("Y", (n, n))

        self.eqns = Equations(
                        Equal(X, Times(A, B, C)),
                        Equal(Y, Times(B, C, D))
                        )


class Example033():
    def __init__(self):

        # TAGS
        # CSE plus

        n = 10

        A = Matrix("A", (n, n))

        B = Matrix("B", (n, n))

        C = Matrix("C", (n, n))

        D = Matrix("D", (n, n))

        X = Matrix("X", (n, n))

        Y = Matrix("Y", (n, n))

        self.eqns = Equations(
                        Equal(X, Plus(A, B, C)),
                        Equal(Y, Plus(B, C, D))
                        )


class Example034():
    def __init__(self):

        # TAGS
        # CSE plus

        n = 10

        A = Matrix("A", (n, n))

        B = Matrix("B", (n, n))

        C = Matrix("C", (n, n))

        D = Matrix("D", (n, n))

        X = Matrix("X", (n, n))

        Y = Matrix("Y", (n, n))

        Z = Matrix("Z", (n, n))
        Z.set_property(properties.AUXILIARY)

        self.eqns = Equations(
                        Equal(Z, Plus(A, B)),
                        Equal(X, Plus(C, Z)),
                        Equal(Y, Plus(Z, D))
                        )


class Example035():
    def __init__(self):

        # TAGS
        # matrix chain, simple

        n1 = 10
        n2 = 20
        n3 = 15
        n4 = 30

        A = Matrix("A", (n1, n2))

        B = Matrix("B", (n2, n3))

        C = Matrix("C", (n3, n4))

        d = Vector("d", (n4, 1))

        x = Vector("x", (n1, 1))

        self.eqns = Equations(
                        Equal(x, Times(A, B, C, d))
                        )
        

class Example036():
    def __init__(self):

        # TAGS
        # random, complex inverse

        n1 = 10
        n2 = 20
        n3 = 15

        A = Matrix("A", (n1, n1))
        A.set_property(properties.SPD)

        B = Matrix("B", (n1, n2))

        C = Matrix("C", (n2, n1))

        D = Matrix("D", (n1, n3))

        X = Matrix("X", (n1, n3))

        self.eqns = Equations(
                        Equal(X, Times(Inverse(Times(A, B, C)), D))
                        )


class Example037():
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

        B = Matrix("B", (n2, n3))

        C = Matrix("C", (n3, n4))

        d = Vector("d", (n4, 1))

        x = Vector("x", (n1, 1))

        E = Matrix("E", (n5, n6))

        F = Matrix("F", (n5, n6))

        G = Matrix("G", (n5, n6))

        Y = Matrix("Y", (n5, n6))

        self.eqns = Equations(
                        Equal(x, Times(A, B, C, d)),
                        Equal(Y, Plus(E, F, G))
                        )


class Example038():
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


        X = Matrix("X", (n1, n7), set([i, j ,k]))


        self.eqns = Equations(
                            # Equal(X, Times(A, B, C, D, E, F, G, H))
                            Equal(X, Times(A, B, C, D, E, F)) 
                            )




class Example039():
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


        A.set_property(properties.LOWER_TRIANGULAR)
        B.set_property(properties.LOWER_TRIANGULAR)
        F.set_property(properties.LOWER_TRIANGULAR)


        X = Matrix("X", (n1, n2), set([i, j ,k]))

        self.eqns = Equations(
                            # Equal(X, Plus(A, B, C, D, E, F, G, H)) 
                            Equal(X, Plus(A, B, C, D, E, F)) 
                            )
        

class Example040():
    def __init__(self):

        # TAGS
        # CSE times

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))


        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Plus(A, B, C, Transpose(A), Transpose(B))),
                            Equal(Y, Plus(A, B))
                            )


class Example041():
    def __init__(self):

        # TAGS
        # CSE plus

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        D = Matrix("D", (n1, n1))


        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        Z = Matrix("Z", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Plus(Transpose(A), B, C,)),
                            Equal(Y, Plus(B, C, D)),
                            Equal(Z, Plus(Transpose(B), A, C))
                            )


class Example042():
    def __init__(self):

        # TAGS
        # GWAS, generalized least squares (inverse missing), complex inverse
        # symmetric product

        n = 10
        m = 5

        S = Matrix("S", (n, n))
        S.set_property(properties.SPD)

        CP = Matrix("CP", (n, m))
        CP.set_property(properties.COLUMN_PANEL)
        CP.set_property(properties.FULL_RANK)

        z = Vector("z", (m, 1))

        y = Vector("y", (n, 1))

        self.eqns = Equations(Equal(z, Times(Inverse(Times(Transpose(CP), S, CP ) ), Transpose(CP), Inverse(S), y)))


class Example043():
    def __init__(self):

        # TAGS
        # CSE times

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))

        A.set_property(properties.LOWER_TRIANGULAR)
        A.set_property(properties.SQUARE)
        B.set_property(properties.LOWER_TRIANGULAR)
        B.set_property(properties.SQUARE)

        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Times(A, B, C, Transpose(Inverse(A)), Transpose(Inverse(B)))),
                            Equal(Y, Times(A, B))
                            )


class Example044():
    def __init__(self):

        # TAGS
        # CSE times

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))

        A.set_property(properties.LOWER_TRIANGULAR)
        A.set_property(properties.SQUARE)
        B.set_property(properties.LOWER_TRIANGULAR)
        B.set_property(properties.SQUARE)

        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Times(Transpose(B), Transpose(A), C, Transpose(Inverse(A)), Transpose(Inverse(B)))),
                            Equal(Y, Times(A, B))
                            )

class Example045():
    def __init__(self):

        # TAGS
        # CSE times

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))


        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Times(A, B, C, Transpose(Inverse(A)), Transpose(Inverse(B)))),
                            Equal(Y, Times(A, B))
                            )


class Example046():
    def __init__(self):

        # TAGS
        # CSE times

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))


        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Times(A, B, C, Transpose(Inverse(A)), Transpose(Inverse(B)))),
                            Equal(Y, Times(A, B))
                            )


class Example047():
    def __init__(self):

        # TAGS
        # POS, SOP

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        D = Matrix("D", (n1, n1))


        X = Matrix("X", (n1, n1))

        self.eqns = Equations(
                            # Equal(X, Times(Plus(A, B), Plus(C, Times(D, Plus(B, C)))))
                            # Equal(X, Times(Plus(A, B), Plus(C, D)))
                            # Equal(X, Times(Plus(A, B), C))
                            Equal(X, Plus(Times(A, Plus(C, D)), Times(B, Plus(C, D))))
                            )


class Example048():
    def __init__(self):

        # TAGS
        # POS, SOP

        n1 = 10
        n2 = 10 

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        D = Matrix("D", (n1, n1))


        X = Matrix("X", (n1, n1))

        # from patternmatcher.expression import to_SOP

        self.eqns = Equations(
                            # Equal(X, Times(Plus(A, B), Plus(C, Times(D, Plus(B, C)))))
                            Equal(X, Times(Plus(A, B), Plus(C, D)))
                            # Equal(X, to_SOP(Times(Plus(A, Transpose(Inverse(Plus(A, B)))), Plus(C, Inverse(Plus(A, B))))))
                            # Equal(X, Times(Plus(A, B), C))
                            # Equal(X, Plus(Times(A, Plus(B, C)), Times(D, Plus(B, C))))
                            )


class Example049():
    def __init__(self):

        # TAGS
        # linear system, matrix chain

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))

        A.set_property(properties.NON_SINGULAR)

        X = Matrix("X", (n1, n1))

        # from patternmatcher.expression import to_SOP

        self.eqns = Equations(
                            Equal(X, Times(InverseTranspose(A), B))
                            )


class Example050():
    def __init__(self):

        # TAGS
        # CSE plus times, random

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))


        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Plus(Times(A, B), C, Inverse(Times(A, B))))
                            )


class Example051():
    def __init__(self):

        # TAGS
        # CSE, random, hard

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))


        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        Z = Matrix("Z", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Plus(A, C, Inverse(Times(A, B)))),
                            Equal(Y, Plus(B, C, Inverse(Times(A, B)))),
                            Equal(Z, Plus(B, C, InverseTranspose(Times(A, B))))
                            )


class Example052():
    def __init__(self):

        # TAGS
        # CSE plus times, random, hard

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))


        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        Z = Matrix("Z", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Plus(A, C, Inverse(Times(A, B)))),
                            Equal(Y, Plus(A, C, Inverse(Times(A, C)))),
                            Equal(Y, Plus(B, C, Inverse(Times(A, C)))),
                            Equal(Z, Plus(B, C, InverseTranspose(Times(A, B))))
                            )


class Example053():
    def __init__(self):

        # TAGS
        # CSE times

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        D = Matrix("D", (n1, n1))


        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        Z = Matrix("Z", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Times(A, B, C, D)),
                            Equal(Y, Times(C, D, A, B))
                            )

class Example054():
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


        X = Matrix("X", (n1, n1))

        self.eqns = Equations(
                            # Equal(X, Times(Plus(A, B), Plus(C, Times(D, Plus(B, C)))))
                            # Equal(X, Times(Plus(A, B), Plus(C, D)))
                            # Equal(X, Times(Plus(A, B), C))
                            Equal(X, Plus(Times(A, Plus(C, D, E)), Times(B, Plus(C, D))))
                            )


class Example055():
    def __init__(self):

        # TAGS
        # CSE times, random

        n1 = 10
        n2 = 10
        n3 = 10

        A = Matrix("A", (n1, n1))
        A.set_property(properties.SPD)

        B = Matrix("B", (n1, n2))

        C = Matrix("C", (n2, n1))

        D = Matrix("D", (n1, n3))

        X = Matrix("X", (n1, n3))

        Y = Matrix("Y", (n1, n1))

        self.eqns = Equations(
                        Equal(X, Times(Inverse(Times(A, B, C)), D)),
                        Equal(Y, Times(B, C))
                        )


class Example056():
    def __init__(self):

        # TAGS
        # CSE plus

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        D = Matrix("D", (n1, n1))


        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        Z = Matrix("Z", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Plus(A, B, C,)),
                            Equal(Y, Plus(B, C, D)),
                            Equal(Z, Plus(D, A, C))
                            )


class Example057():
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


        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        Z = Matrix("Z", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Plus(A, B, C, Transpose(A), Transpose(B))),
                            # Equal(X, Plus(A, B, D, Transpose(A), Transpose(B), E, C)),
                            # Equal(Y, Plus(A, B, D, Transpose(A), Transpose(B))),
                            # Equal(Y, Plus(A, B, D, Transpose(A), Transpose(B), E)),
                            # Equal(Z, Plus(A, B)),
                            # Equal(Z, Plus(D, E)),
                            # Equal(X, Plus(C, Transpose(A), Transpose(B)))
                            )


class Example058():
    def __init__(self):

        # TAGS
        # CSE plus

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        D = Matrix("D", (n1, n1))

        D.set_property(properties.DIAGONAL)

        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        Z = Matrix("Z", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Plus(A, B, D)),
                            Equal(Y, Plus(C, B, D)),
                            Equal(Z, Plus(A, B)),
                            )


class Example059():
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


        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        Z = Matrix("Z", (n1, n1))

        # X = (A+D)B(C+E)
        self.eqns = Equations(
                            Equal(X, Times(Plus(A, D), B, Plus(C, E))),
                            # Equal(X, Plus(Times(A, B, C), Times(D, B, C), Times(D, B, E), Times(A, B, E))),
                            )


class Example060():
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


        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        Z = Matrix("Z", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Plus(Times(A, B), Times(A, C))),
                            )


class Example061():
    def __init__(self):

        # TAGS
        # least squares, QR

        n = 10
        m = 5

        X = Matrix("X", (n, m))
        X.set_property(properties.FULL_RANK)
        X.set_property(properties.COLUMN_PANEL)

        y = Vector("y", (n, 1))

        b = Vector("b", (m, 1))

        # b = (X^T X)^-1 X^T y
        self.eqns = Equations(Equal(b, Times(Inverse(Times(Transpose(X), X)), Transpose(X), y)))


class Example062():
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


        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        Z = Matrix("Z", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Plus(Times(alpha, A), Times(beta, B), C)),
                            )


class Example063():
    def __init__(self):

        # TAGS
        # GWAS, generalized least squares

        n = 10 # 10
        m = 5 # 5

        S = Matrix("S", (n, n))
        S.set_property(properties.SPD)

        X = Matrix("X", (n, m))
        # X.set_property(properties.COLUMN_PANEL)
        X.set_property(properties.FULL_RANK)

        z = Vector("z", (m, 1))

        y = Vector("y", (n, 1))

        # clak.special_properties.add_expression(Times(Transpose(X), Inverse(S), X ), set([properties.SPD))

        self.eqns = Equations(Equal(z, Times(Inverse(Times(Transpose(X), Inverse(S), X ) ), Transpose(X), Inverse(S), y)))


class Example064():
    def __init__(self):

        # TAGS
        # transposed kernels

        n = 10

        A = Matrix("A", (n, n))

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)

        x = Vector("x", (n, 1))

        y = Vector("y", (n, 1))

        w = Vector("w", (1, n))

        self.eqns = Equations(
                        Equal(w, Times(Transpose(x), B)),
                        )


class Example065():
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


        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        Z = Matrix("Z", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Times(A, B, C, D)),
                            Equal(Y, Times(B, C, D, E)),
                            )


class Example066():
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


        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        Z = Matrix("Z", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Times(A, B, C, D)),
                            Equal(Y, Times(B, C, D, E)),
                            Equal(Z, Times(B, C, D, F)),
                            )


class Example067():
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


        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        Z = Matrix("Z", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Times(A, B, C, D)),
                            Equal(Y, Times(B, C, D, E)),
                            Equal(Z, Times(B, C, D, E, F)),
                            )


class Example068():
    def __init__(self):

        # TAGS
        # matrix chain

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))

        A.set_property(properties.LOWER_TRIANGULAR)
        A.set_property(properties.SQUARE)
        B.set_property(properties.LOWER_TRIANGULAR)
        B.set_property(properties.SQUARE)

        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        self.eqns = Equations(
                            Equal(Y, Times(C, InverseTranspose(A), InverseTranspose(B)))
                            )


class Example069():
    def __init__(self):

        # TAGS
        # CSE times

        n1 = 10
        n2 = 20

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))

        A.set_property(properties.LOWER_TRIANGULAR)
        A.set_property(properties.SQUARE)
        B.set_property(properties.LOWER_TRIANGULAR)
        B.set_property(properties.SQUARE)

        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Times(Transpose(B), Transpose(A), C, Transpose(Inverse(A)), Transpose(Inverse(B)))),
                            )


class Example070():
    def __init__(self):

        # TAGS
        # complex inverse, random, explicit inverse

        n = 10

        A = Matrix("A", (n, n))

        D = Matrix("D", (n, n))
        D.set_property(properties.DIAGONAL)

        W = Matrix("W", (n, n))

        self.eqns = Equations(Equal(W, Times(Inverse(Plus(A, Inverse(D) ) ), D)))


class Example071():
    def __init__(self):

        # TAGS
        # complex inverse, random

        n = 10

        A = Matrix("A", (n, n))

        D = Matrix("D", (n, n))
        D.set_property(properties.DIAGONAL)

        W = Matrix("W", (n, n))

        self.eqns = Equations(Equal(W, Times(Inverse(Plus(A, D) ), D)))


class Example072():
    def __init__(self):

        # TAGS
        # complex inverse, termination/progress is a problem

        n = 10

        A = Matrix("A", (n, n))

        B = Matrix("B", (n, n))

        W = Matrix("W", (n, n))

        self.eqns = Equations(Equal(W, Plus(Inverse(A), B, A)))


class Example073():
    def __init__(self):

        # TAGS
        # explicit inverse, simple

        n = 10

        A = Matrix("A", (n, n))

        D = Matrix("D", (n, n))
        D.set_property(properties.DIAGONAL)

        W = Matrix("W", (n, n))

        self.eqns = Equations(Equal(W, Plus(Inverse(D), A)))


class Example074():
    def __init__(self):

        # TAGS
        # eigen trick

        n = 10

        Q = Matrix("Q", (n, n))
        Q.set_property(properties.ORTHOGONAL)

        W = Matrix("W", (n, n))
        W.set_property(properties.DIAGONAL)

        I = IdentityMatrix(n, n)

        alpha = Scalar("alpha")

        two = ConstantScalar(2.0)

        X = Matrix("X", (n, n))

        self.eqns = Equations(Equal(X, Plus(Times(Transpose(Q), W, Q), Times(Plus(two, alpha), I))))


class Example075():
    def __init__(self):

        # TAGS
        # canonical form, simplify, SOP, POS, scalar in sum

        n = 10

        A = Matrix("A", (n, n))

        alpha = Scalar("alpha")

        two = ConstantScalar(2.0)

        X = Matrix("X", (n, n))

        self.eqns = Equations(Equal(X, Plus(Times(alpha, A), Times(two, A))))


class Example076():
    def __init__(self):

        # TAGS
        # eigen trick, SOP, POS

        n = 10

        A = Matrix("A", (n, n))

        B = Matrix("B", (n, n))

        Q = Matrix("Q", (n, n))
        Q.set_property(properties.ORTHOGONAL)

        W = Matrix("W", (n, n))
        W.set_property(properties.DIAGONAL)

        I = IdentityMatrix(n, n)

        alpha = Scalar("alpha")

        two = ConstantScalar(2.0)

        X = Matrix("X", (n, n))

        self.eqns = Equations(Equal(X, Times(A, Plus(Times(Transpose(Q), W, Q), Times(alpha, I)), B)))


class Example077():
    def __init__(self):

        # TAGS
        # eigen trick, equivalent expression, SOP, POS

        n = 10

        A = Matrix("A", (n, n))

        B = Matrix("B", (n, n))

        Q = Matrix("Q", (n, n))
        Q.set_property(properties.ORTHOGONAL)

        W = Matrix("W", (n, n))
        W.set_property(properties.DIAGONAL)

        I = IdentityMatrix(n, n)

        alpha = Scalar("alpha")

        two = ConstantScalar(2.0)

        X = Matrix("X", (n, n))

        Y = Matrix("Y", (n, n))

        self.eqns = Equations(
                        Equal(Y, Plus(A, B)),
                        Equal(X, Times(A, Plus(Times(Transpose(Q), W, Q), Times(alpha, I)), B))
                        )


class Example078():
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


        S.set_property(properties.SPD)

        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        Z = Matrix("Z", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Times(A, B, S, Transpose(B), Transpose(A))),
                            )


class Example079():
    def __init__(self):

        # TAGS
        # CSE times

        n1 = 10
        n2 = 20

        A = Matrix("A", (n2, n1))
        B = Matrix("B", (n1, n1))


        X = Matrix("X", (n2, n1))

        self.eqns = Equations(
                            Equal(X, Times(Inverse(Times(A, B, Transpose(B), Transpose(A))), A, B)),
                            )


class Example080():
    def __init__(self):

        # TAGS
        # CSE times

        n1 = 10
        n2 = 20

        A = Matrix("A", (n2, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n2, n2))


        X = Matrix("X", (n2, n1))

        Y = Matrix("Y", (n2, n2))

        self.eqns = Equations(
                            Equal(X, Times(C, A, B)),
                            Equal(Y, Times(A, B, Transpose(B), Transpose(A)))
                            )


class Example081():
    def __init__(self):

        # TAGS
        # CSE times

        n = 10

        A = Matrix("A", (n, n))

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)

        C = Matrix("C", (n, n))

        W = Matrix("W", (n, n))

        self.eqns = Equations(Equal(W, Times(InverseTranspose(B), InverseTranspose(C), A, Transpose(C), Transpose(B))))


class Example082():
    def __init__(self):

        # TAGS
        # CSE times

        n = 10

        A = Matrix("A", (n, n))

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)

        C = Matrix("C", (n, n))

        W = Matrix("W", (n, n))

        self.eqns = Equations(Equal(W, Times(Transpose(B), InverseTranspose(C), A, Transpose(C), InverseTranspose(B))))


class Example083():
    def __init__(self):

        # TAGS
        # CSE times

        n = 10

        A = Matrix("A", (n, n))
        D = Matrix("D", (n, n))

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)

        C = Matrix("C", (n, n))

        W = Matrix("W", (n, n))

        self.eqns = Equations(Equal(W, Times(InverseTranspose(B), InverseTranspose(C), A, Transpose(C), Transpose(B),  D, Inverse(C), Inverse(B))))


class Example084():
    def __init__(self):

        # TAGS
        # CSE times, termination/progress is a problem

        n = 10

        A = Matrix("A", (n, n))

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)

        C = Matrix("C", (n, n))

        W = Matrix("W", (n, n))

        self.eqns = Equations(Equal(W, Plus(Times(B, C), Times(Inverse(B), A))))


class Example085():
    def __init__(self):

        # TAGS
        # CSE times, symmetric product

        n1 = 10
        n2 = 10

        A = Matrix("A", (n2, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n2, n2))


        X = Matrix("X", (n2, n1))

        Y = Matrix("Y", (n1, n1))

        self.eqns = Equations(
                            Equal(Y, Times(Transpose(A), Transpose(B), B, A))
                            )

class Example086():
    def __init__(self):

        # TAGS
        # CSE times, symmetric product

        n1 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))


        A.set_property(properties.LOWER_TRIANGULAR)
        B.set_property(properties.LOWER_TRIANGULAR)
        C.set_property(properties.LOWER_TRIANGULAR)

        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Plus(Times(A, B, Transpose(B), Transpose(A)), C))
                            )

class Example087():
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


        X = Matrix("X", (n1, n5), set([i, j]))


        self.eqns = Equations(
                            Equal(X, Times(A, B, C, D)) 
                            )

class Example088():
    def __init__(self):

        # TAGS
        # trick, tricks

        n1 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))


        X = Matrix("X", (n1, n1))


        self.eqns = Equations(
                            Equal(X, Plus(Times(Transpose(A), B), Times(Transpose(B), A), Times(Transpose(A), A), C)) 
                            )

class Example089():
    def __init__(self):

        # TAGS
        # eigen trick

        n = 10

        A = Matrix("A", (n, n))

        B = Matrix("B", (n, n))

        Q = Matrix("Q", (n, n))
        Q.set_property(properties.ORTHOGONAL)

        W = Matrix("W", (n, n))
        W.set_property(properties.DIAGONAL)

        I = IdentityMatrix(n, n)

        alpha = Scalar("alpha")

        two = ConstantScalar(2.0)

        X = Matrix("X", (n, n))

        self.eqns = Equations(Equal(X, Plus(Times(Transpose(Q), W, Q), Times(alpha, I), B)))


class Example090():
    def __init__(self):

        # TAGS
        # blocks, blocked, partitioning

        n = 10

        A = Matrix("A", (n, n))

        B = Matrix("B", (n, n))
        B.partitioning = (set([5]), set([5]))

        D = Matrix("D", (n, n))
        # D.set_property(properties.DIAGONAL)

        alpha = Scalar("alpha")

        I = IdentityMatrix(n, n)

        X = Matrix("X", (n, n))

        self.eqns = Equations(Equal(X, Times(alpha, B, D)))


class Example091():
    def __init__(self):

        # TAGS
        # blocks, blocked, partitioning, indices

        n = 10

        i = Index("i", 50)

        A = Matrix("A", (n, n))

        B = Matrix("B", (n, n), set([i]))

        C = Matrix("C", (2*n, n))
        # D.set_property(properties.DIAGONAL)

        alpha = Scalar("alpha")

        I = IdentityMatrix(n, n)

        X = Matrix("X", (n, n), set([i]))

        self.eqns = Equations(Equal(X, Times(BlockedExpression([[A, B]]), C)))


class Example092():
    def __init__(self):

        # TAGS
        # explicit inverse, matrix chain

        n = 100
        m = 10

        V = Matrix("V", (n, m))
        V.set_property(properties.FULL_RANK)

        X = Matrix("X", (m, m))


        self.eqns = Equations(Equal(X, Inverse(Times(Transpose(V), V))))
        # self.eqns = Equations(Equal(X, Times(Transpose(V), V)))


class Example093():
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


        self.eqns = Equations(Equal(X, Times(A, Inverse(B), Transpose(A), C)),
                              Equal(X, Times(D, A, Inverse(B), Transpose(A)))
                              )


class Example094():
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

        A.set_property(properties.SPD)

        X = Matrix("X", (n1, n5), set([i, j]))


        self.eqns = Equations(
                            # Equal(X, Times(A, B, C, D, E, F, G, H))
                            Equal(X, Times(Inverse(A), B, C, D)) 
                            )

class Example095():
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



        X = Matrix("X", (n1, 1), set([i, j]))


        self.eqns = Equations(
                            Equal(X, Times(A, B, c)) 
                            )

class Example096():
    def __init__(self):

        # TAGS
        # simple, matrix chain

        n1 = 1500
        n2 = 1000
        n3 = 500

        A = Matrix("A", (n1, n2))

        B = Matrix("B", (n2, n2))
        B.set_property(properties.LOWER_TRIANGULAR)

        C = Matrix("C", (n2, n2))
        C.set_property(properties.UPPER_TRIANGULAR)

        D = Matrix("D", (n2, n3))
        D.set_property(properties.DIAGONAL)

        W = Matrix("W", (n1, n3))

        self.eqns = Equations(Equal(W, Times(A, Inverse(B), InverseTranspose(C), D)))


class Example097():
    def __init__(self):

        # TAGS
        # simple, matrix chain

        n1 = 600
        n2 = 200
        n3 = 1500
        n4 = 800
        n5 = 1000

        # n1 = 1500
        # n2 = 1200
        # n3 = 1300
        # n4 = 800
        # n5 = 1000

        A = Matrix("A", (n1, n1))
        A.set_property(properties.SPD)

        B = Matrix("B", (n1, n2))

        C = Matrix("C", (n2, n3))

        D = Matrix("D", (n3, n3))
        D.set_property(properties.UPPER_TRIANGULAR)

        E = Matrix("E", (n3, n4))
        E.set_property(properties.UPPER_TRIANGULAR)

        F = Matrix("F", (n4, n5))

        W = Matrix("W", (n1, n5))

        self.eqns = Equations(Equal(W, Times(Inverse(A), B, C, InverseTranspose(D), E, F)))


class Example098():
    def __init__(self):

        # TAGS
        # simple, scalar, constant

        alpha = Scalar("alpha")
        two = ConstantScalar(2.0)

        beta = Scalar("beta")

        self.eqns = Equations(Equal(beta, Times(alpha, two)))


class Example099():
    def __init__(self):

        # TAGS
        # simple, scalar, constant

        alpha = Scalar("alpha")
        two = ConstantScalar(2.0)

        beta = Scalar("beta")

        self.eqns = Equations(Equal(beta, Plus(alpha, two)))


class Example100():
    def __init__(self):

        # TAGS
        # simple, scalar, constant

        n = 10

        A = Matrix("A", (n, n))
        two = ConstantScalar(2.0)

        X = Matrix("X", (n, n))


        beta = Scalar("beta")

        self.eqns = Equations(Equal(X, Times(two, A)))


class Example101():
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


class Example102():
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


class Example103():
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


class Example104():
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



class Example105():
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


class Example106():
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


class Example107():
    def __init__(self):

        # TAGS
        # simple, plus, identity

        n1 = 50

        A = Matrix("A", (n1, n1))
        K = Matrix("K", (n1, n1))
        X = Matrix("X", (n1, n1))


        self.eqns = Equations(Equal(X, Plus(A, IdentityMatrix(n1, n1))))


class Example108():
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


class Example109():
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


class Example110():
    def __init__(self):

        # TAGS
        # eigen trick

        n = 10

        Q = Matrix("Q", (n, n))
        Q.set_property(properties.ORTHOGONAL)

        W = Matrix("W", (n, n))
        W.set_property(properties.DIAGONAL)

        I = IdentityMatrix(n, n)

        alpha = Scalar("alpha")

        two = ConstantScalar(2.0)

        X = Matrix("X", (n, n))

        self.eqns = Equations(Equal(X, Plus(Times(Transpose(Q), W, Q), Times(alpha, I))))


class Example111():
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


class Example112():
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


class Example113():
    def __init__(self):

        # TAGS
        # POS, SOP

        n1 = 10
        n2 = 10 

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        D = Matrix("D", (n1, n1))


        X = Matrix("X", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Times(Plus(A, B), C))
                            )


class Example114():
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


class Example115():
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

        B = Matrix("B", (n1, n2))

        C = Matrix("C", (n2, n3))

        D = Matrix("D", (n3, n3))
        D.set_property(properties.UPPER_TRIANGULAR)

        E = Matrix("E", (n3, n4))
        E.set_property(properties.UPPER_TRIANGULAR)

        F = Matrix("F", (n4, n5))

        W = Matrix("W", (n1, n5))

        self.eqns = Equations(Equal(W, Times(A, B, C, D, E, F)))


class Example116():
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


        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        Z = Matrix("Z", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Times(A, Transpose(B), C, B, Transpose(A))),
                            )


class Example117():
    def __init__(self):

        # TAGS
        # matrix chain, linear system

        n = 10

        S = Matrix("B", (n, n))
        S.set_property(properties.SPD)

        x = Vector("x", (n, 1))

        w = Vector("w", (n, 1))

        self.eqns = Equations(Equal(w, Times(Inverse(S), x)))


class Example118():
    def __init__(self):

        # TAGS
        # simple, sum, inverse

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.NON_SINGULAR)

        X = Matrix("X", (n, n))

        self.eqns = Equations(Equal(X, Plus(A, Inverse(A))))


class Example119():
    def __init__(self):

        # TAGS
        # simple, sum, inverse

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.SPD)

        X = Matrix("X", (n, n))

        self.eqns = Equations(Equal(X, Plus(A, Inverse(A))))


class Example120():
    def __init__(self):

        # TAGS
        # trick, tricks

        n = 20
        m = 10

        A = Matrix("A", (n, m))
        B = Matrix("B", (n, m))


        X = Matrix("X", (m, m))


        self.eqns = Equations(
                            Equal(X, Plus(Times(Transpose(A), B), Times(Transpose(B), A), Times(Transpose(A), A))) 
                            )


class Example121():
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


class Example122():
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


class Example123():
    def __init__(self):

        # TAGS
        # noschese2016

        n = 1500 # 10
        m = 1000 # 5

        A = Matrix("A", (n, m))
        A.set_property(properties.FULL_RANK)

        b = Vector("b", (n, 1))

        mu = Scalar("mu")

        I = IdentityMatrix(m, m)

        x = Vector("x", (m, 1))

        self.eqns = Equations(
                            Equal(x, 
                                Times(Inverse(Plus(Times(Transpose(A), A), Times(mu, I))), Transpose(A), b)
                                ) 
                            )

class Example124():
    def __init__(self):

        # TAGS
        # CSE times

        n1 = 1000
        n2 = 1000

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        D = Matrix("D", (n1, n1))

        B.set_property(properties.NON_SINGULAR)
        # B.set_property(properties.SPD)
        B.set_property(properties.LOWER_TRIANGULAR)

        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Times(A, InverseTranspose(B), C)),
                            Equal(Y, Times(D, Inverse(B), Transpose(A))),
                            )


class Example125():
    def __init__(self):

        # TAGS
        # nino2016, application

        nsd = 1000
        msd = 1000
        N = 1000

        
        minu1 = ConstantScalar(-1.0)

        B = Matrix("B", (N, N))
        H = Matrix("H", (msd, N))
        R = Matrix("R", (msd, msd))
        Y = Matrix("Y", (msd, N))
        Xb = Matrix("Xb", (nsd, N)) # originally X^b
        X = Matrix("X", (nsd, N))


        self.eqns = Equations(
                            Equal(
                                X,
                                Plus(
                                    Xb,
                                    Times(
                                        Inverse(Plus(Inverse(B), Times(Transpose(H), Inverse(R), H))),
                                        Plus(Y, Times(H, Xb, minu1))
                                        )
                                    )
                                )
                            )

        # self.eqns = Equations(
        #                     Equal(
        #                         X,
        #                         Times(
        #                                 Inverse(Plus(Inverse(B), Times(Transpose(H), Inverse(R), H))),
        #                                 Plus(Y, Times(H, Xb))
        #                                 )
        #                             )
        #                     )
        # print(self.eqns)


class Example126():
    def __init__(self):

        # TAGS
        # straszak2015, application

        # n > m
        m = 1000
        n = 1500

        A = Matrix("A", (m, n))
        A.set_property(properties.FULL_RANK)

        W = Matrix("W", (n, n))
        # C = Matrix("C", (n1, n1))

        # W is positive
        W.set_property(properties.FULL_RANK)
        W.set_property(properties.DIAGONAL)
        W.set_property(properties.SPD)
        W.set_property(properties.SYMMETRIC)
        W.set_property(properties.NON_SINGULAR)

        b = Vector("b", (m, 1))

        c = Vector("c", (n, 1))

        x = Vector("x", (n, 1))

        minu1 = ConstantScalar(-1.0)

        # derivation.special_properties.add_expression(Times(A, W, Transpose(A)), {properties.SPD})

        self.eqns = Equations(
                            Equal(
                                x,
                                Times(
                                    W,
                                    Plus(
                                        Times(Transpose(A), Inverse(Times(A, W, Transpose(A))), b),
                                        Times(minu1, c)
                                        )
                                    )
                                )
                            )

class Example127():
    def __init__(self):

        # TAGS
        # trick, tricks

        n1 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        D = Matrix("D", (n1, n1))

        D.set_property(properties.SYMMETRIC)

        X = Matrix("X", (n1, n1))


        self.eqns = Equations(
                            Equal(X, Plus(Times(Transpose(A), B), Times(Transpose(B), A), Times(Transpose(A), D, A), C)) 
                            )


class Example128():
    def __init__(self):

        # TAGS

        # n1 = 55
        # n2 = 58
        # n3 = 64
        # n4 = 54
        # n5 = 128
        # n6 = 53

        n1 = 130
        n2 = 700
        n3 = 383
        n4 = 1340
        n5 = 193
        n6 = 900


        A = Matrix("A", (n1, n2))
        B = Matrix("B", (n2, n3))
        C = Matrix("C", (n3, n4))
        D = Matrix("D", (n4, n5))
        E = Matrix("E", (n5, n6))


        X = Matrix("X", (n1, n6))

        self.eqns = Equations(
                            Equal(X, Times(A, B, C, D, E)),
                            )

class Example129():
    def __init__(self):

        # TAGS
        # CSE, explicit inversion

        n = 10

        A = Matrix("A", (n, n))
        B = Matrix("B", (n, n))
        C = Matrix("C", (n, n))

        A.set_property(properties.NON_SINGULAR)

        X = Matrix("X", (n, n))

        Y = Matrix("Y", (n, n))

        self.eqns = Equations(
                            Equal(X, Plus(Inverse(A), B)),
                            Equal(Y, Plus(Inverse(A), C)),
                            )

class Example130():
    def __init__(self):

        # TAGS
        # CSE, explicit inversion

        n = 10

        A = Matrix("A", (n, n))
        B = Matrix("B", (n, n))
        C = Matrix("C", (n, n))

        A.set_property(properties.NON_SINGULAR)

        X = Matrix("X", (n, n))

        self.eqns = Equations(
                            Equal(X, Plus(Times(Inverse(A), B), A, Inverse(A))),
                            )


class Example131():
    def __init__(self):

        # TAGS
        # CSE, explicit inversion

        n = 10

        A = Matrix("A", (n, n))
        B = Matrix("B", (n, n))
        C = Matrix("C", (n, n))
        D = Matrix("D", (n, n))

        A.set_property(properties.NON_SINGULAR)

        X = Matrix("X", (n, n))

        self.eqns = Equations(
                            Equal(X, Plus(Times(Inverse(A), B), A, Times(A, C), Inverse(A))),
                            )


class Example132():
    def __init__(self):

        # TAGS
        # CSE, explicit inversion

        n = 10

        A = Matrix("A", (n, n))
        B = Matrix("B", (n, n))
        C = Matrix("C", (n, n))
        D = Matrix("D", (n, n))

        A.set_property(properties.NON_SINGULAR)

        X = Matrix("X", (n, n))

        self.eqns = Equations(
                            Equal(X, Plus(Times(D, Inverse(A), B), A, Times(A, C), Inverse(A))),
                            )


class Example133():
    def __init__(self):

        # TAGS
        # undistribute inverse

        n = 10

        A = Matrix("A", (n, n))
        B = Matrix("B", (n, n))
        C = Matrix("C", (n, n))
        D = Matrix("D", (n, n))

        A.set_property(properties.NON_SINGULAR)
        B.set_property(properties.NON_SINGULAR)

        X = Matrix("X", (n, n))

        self.eqns = Equations(
                            Equal(X, Times(Inverse(A), InverseTranspose(B), C)),
                            )


class Example134():
    def __init__(self):

        # TAGS
        # least squares, CSE

        n = 1500
        m = 1000

        A = Matrix("A", (m, m))
        A.set_property(properties.FULL_RANK)

        B = Matrix("B", (m, m))
        B.set_property(properties.FULL_RANK)

        X = Matrix("X", (n, m))
        X.set_property(properties.FULL_RANK)

        y = Vector("y", (n, 1))

        z = Vector("z", (n, 1))

        b = Vector("b", (m, 1))

        self.eqns = Equations(Equal(b, 
                        Plus(
                            Times(A, Inverse(Times(Transpose(X), X)), Transpose(X), y),
                            Times(B, Inverse(Times(Transpose(X), X)), Transpose(X), z),
                            )
                        ))


class Example135():
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


        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        Z = Matrix("Z", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Times(A, B, C, D)),
                            Equal(Y, Times(A, B, C, E)),
                            Equal(Z, Times(F, B, C, D)),
                            )


class Example136():
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


        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        Z = Matrix("Z", (n1, n1))

        W = Matrix("W", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Times(A, B, C, D)),
                            Equal(Y, Times(A, B, C, E)),
                            Equal(Z, Times(F, B, C, D)),
                            Equal(W, Times(A, B, C, D)),
                            )

class Example137():
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


        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        Z = Matrix("Z", (n1, n1))

        W = Matrix("W", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Times(A, B, C)),
                            Equal(Y, Times(A, B, C)),
                            Equal(Z, Times(B, C, F)),
                            Equal(W, Times(A, B, C)),
                            )

class Example138():
    def __init__(self):

        # TAGS
        # CSE times

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))


        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        Z = Matrix("Z", (n1, n1))

        W = Matrix("W", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Times(A, B, C)),
                            Equal(Y, Times(A, B, C)),
                            Equal(Z, Times(A, B)),
                            Equal(W, Times(A, B)),
                            )

class Example139():
    def __init__(self):

        # TAGS
        # CSE times

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        D = Matrix("D", (n1, n1))


        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        Z = Matrix("Z", (n1, n1))

        W = Matrix("W", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Times(A, B, C, D)),
                            Equal(Y, Times(A, B, C, D)),
                            Equal(Z, Times(A, B)),
                            Equal(W, Times(A, B)),
                            )

class Example140():
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


        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        Z = Matrix("Z", (n1, n1))

        W = Matrix("W", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Times(B, C, D)),
                            Equal(Y, Times(A, B, C, D)),
                            Equal(Z, Times(A, B)),
                            )

class Example141():
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


        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        Z = Matrix("Z", (n1, n1))

        W = Matrix("W", (n1, n1))

        U = Matrix("U", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Times(A, B, C, D, E)),
                            Equal(Y, Times(A, B, C)),
                            Equal(Z, Times(A, B)),
                            Equal(W, Times(C, D, E)),
                            Equal(U, Times(C, D)),
                            )

class Example142():
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
        S = Matrix("S", (n1, n1))

        S.set_property(properties.SYMMETRIC)


        X = Matrix("X", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        Z = Matrix("Z", (n1, n1))

        W = Matrix("W", (n1, n1))

        U = Matrix("U", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Times(A, S, Transpose(A))),
                            Equal(Y, Times(A, S)),
                            Equal(Z, Times(A, B)),
                            Equal(W, Times(A, B, C)),
                            )


class Example143():
    def __init__(self):

        # TAGS
        # linear system, matrix chain

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.SQUARE)
        A.set_property(properties.NON_SINGULAR)

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)

        X = Matrix("X", (n, n))

        self.eqns = Equations(Equal(X, Times(Inverse(A), B)))


class Example144():
    def __init__(self):

        # TAGS
        # slides

        n = 10 # 10
        m = 5 # 5

        M = Matrix("M", (n, n))
        M.set_property(properties.SPD)

        X = Matrix("X", (n, m))
        # X.set_property(properties.COLUMN_PANEL)
        X.set_property(properties.FULL_RANK)

        z = Vector("z", (m, 1))

        y = Vector("y", (m, 1))

        self.eqns = Equations(Equal(z, Times(Transpose(X), Inverse(M), X, y)))


class Example145():
    def __init__(self):

        # TAGS
        # slides

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        D = Matrix("D", (n1, n1))


        X = Matrix("X", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Plus(Times(A, B, C), D)),
                            )


class Example146():
    def __init__(self):

        # TAGS
        # slides

        n1 = 10
        n2 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        D = Matrix("D", (n1, n1))
        E = Matrix("E", (n1, n1))


        X = Matrix("X", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Plus(Times(A, B, C), Times(C, D, E))),
                            )


class Example147():
    def __init__(self):

        # TAGS
        # properties,

        n = 100

        L1 = Matrix("L1", (n, n))
        L1.set_property(properties.LOWER_TRIANGULAR)
        L1.set_property(properties.FULL_RANK)

        L2 = Matrix("L2", (n, n))
        L2.set_property(properties.LOWER_TRIANGULAR)
        L2.set_property(properties.FULL_RANK)

        U = Matrix("U", (n, n))
        U.set_property(properties.UPPER_TRIANGULAR)
        U.set_property(properties.FULL_RANK)

        B = Matrix("B", (n, n))

        X = Matrix("X", (n, n))

        self.eqns = Equations(
                            Equal(X, Times(Inverse(Plus(L1, Times(InverseTranspose(U), L2))), B)),
                            )


class Example148():
    def __init__(self):

        # TAGS
        # properties,

        n = 100

        L = Matrix("L", (n, n))
        L.set_property(properties.LOWER_TRIANGULAR)
        L.set_property(properties.FULL_RANK)

        U1 = Matrix("U1", (n, n))
        U1.set_property(properties.UPPER_TRIANGULAR)
        U1.set_property(properties.FULL_RANK)

        U2 = Matrix("U2", (n, n))
        U2.set_property(properties.UPPER_TRIANGULAR)
        U2.set_property(properties.FULL_RANK)

        B = Matrix("B", (n, n))

        X = Matrix("X", (n, n))

        self.eqns = Equations(
                            Equal(X, Times(Inverse(Plus(L, Times(InverseTranspose(U1), U2))), B)),
                            )


class Example149():
    def __init__(self):

        # TAGS
        # trick, tricks

        n1 = 10

        A = Matrix("A", (n1, n1))
        B = Matrix("B", (n1, n1))
        C = Matrix("C", (n1, n1))
        D = Matrix("D", (n1, n1))


        X = Matrix("X", (n1, n1))


        self.eqns = Equations(
                            Equal(X, Plus(Times(Transpose(A), C, B), Times(Transpose(B), Transpose(C), A), Times(Transpose(A), C, Transpose(C), A))) 
                            )


class Example150():
    def __init__(self):

        # TAGS
        # diagonal

        n = 1000

        D = Matrix("D", (n, n))
        D.set_property(properties.SQUARE)
        D.set_property(properties.DIAGONAL)
        D.set_property(properties.FULL_RANK)

        B = Matrix("B", (n, n))
        B.set_property(properties.SQUARE)
        B.set_property(properties.NON_SINGULAR)

        X = Matrix("X", (n, n))

        self.eqns = Equations(Equal(X, Times(D, B)))


class Example151():
    def __init__(self):

        # TAGS
        # explicit inversion

        n = 10

        A = Matrix("A", (n, n))
        B = Matrix("B", (n, n))
        C = Matrix("C", (n, n))

        A.set_property(properties.NON_SINGULAR)

        X = Matrix("X", (n, n))

        self.eqns = Equations(
                            Equal(X, Plus(A, Inverse(A))),
                            )


class Example152():
    def __init__(self):

        # TAGS
        # simple, product

        n1 = 10
        n2 = 10 

        S = Matrix("S", (n1, n1))
        S.set_property(properties.SYMMETRIC)

        B = Matrix("B", (n1, n1))

        X = Matrix("X", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Times(S, B))
                            )


class Example153():
    def __init__(self):

        # TAGS
        # simple, product

        n1 = 10
        n2 = 10 

        A = Matrix("A", (n1, n1))

        X = Matrix("X", (n1, n1))

        self.eqns = Equations(
                            Equal(X, Times(Transpose(A), A))
                            )


class Example154():
    def __init__(self):

        # TAGS
        # linear system, matrix chain

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.NON_SINGULAR)
        A.set_property(properties.SYMMETRIC)

        B = Matrix("B", (n, n))
        B.set_property(properties.NON_SINGULAR)
        B.set_property(properties.SYMMETRIC)

        X = Matrix("X", (n, n))

        self.eqns = Equations(
                            Equal(X, Times(Inverse(A), B)),
                            # Equal(X, Times(A, Inverse(B)))
                            )


class Example155():
    def __init__(self):

        # TAGS
        # matrix chain

        # M: 1150 by 1100
        # v1: 1150 by 1
        # v2: 1150 by 1
        # v3: 1650 by 1
        # v4: 1650 by 1

        # M’*v*v3’*v4

        n = 10

        M = Matrix("M", (1100, 1150))

        v1 = Vector("v1", (1150, 1))

        v2 = Vector("v2", (1650, 1))

        v3 = Vector("v3", (1650, 1))

        x = Vector("x", (1100, 1))

        self.eqns = Equations(Equal(x, Times(M, v1, Transpose(v2), v3)))


class Example156():
    def __init__(self):

        # TAGS
        # linear system, matrix chain

        n = 1000

        A = Matrix("A", (n, n))
        A.set_property(properties.SPD)

        B = Matrix("B", (n, n))
        B.set_property(properties.FULL_RANK)

        C = Matrix("C", (n, n))
        C.set_property(properties.LOWER_TRIANGULAR)
        C.set_property(properties.FULL_RANK)

        X = Matrix("X", (n, n))

        self.eqns = Equations(
                            Equal(X, Times(Inverse(A), B, Transpose(C))),
                            )


class Example157():
    def __init__(self):

        # TAGS
        # simple, diagonal, sum

        n = 10

        A = Matrix("A", (n, n))
        A.set_property(properties.DIAGONAL)

        B = Matrix("B", (n, n))
        B.set_property(properties.DIAGONAL)

        alpha = Scalar("alpha")

        I = IdentityMatrix(n, n)

        X = Matrix("X", (n, n))


        self.eqns = Equations(Equal(X, Plus(Times(alpha, I), B)))


class Example158():
    def __init__(self):

        # TAGS
        # complex inverse

        n1 = 100
        n2 = 100
        # n2 = 50

        A = Matrix("A", (n1, n2))
        A.set_property(properties.FULL_RANK)

        B = Matrix("B", (n2, n1))
        B.set_property(properties.FULL_RANK)

        C = Matrix("C", (n1, n1))

        X = Matrix("X", (n1, n1))

        self.eqns = Equations(
                        Equal(X, Times(Inverse(Times(A, B)), C))
                        # Equal(X, Inverse(Times(A, B)))
                        )


class Example159():
    def __init__(self):

        # TAGS
        # in/out operands

        n1 = 100
        n2 = 100
        # n2 = 50

        A = Matrix("A", (n1, n2))
        A.set_property(properties.FULL_RANK)

        B = Matrix("B", (n2, n1))
        B.set_property(properties.FULL_RANK)

        C = Matrix("C", (n1, n1))

        Y = Matrix("Y", (n1, n1))

        X = Matrix("X", (n1, n1))

        self.eqns = Equations(
                        Equal(Y, Times(Transpose(A), A)),
                        Equal(X, Times(Inverse(Y), C))
                        )


class Example160():
    def __init__(self):

        # TAGS
        # SAC

        n = 5000

        A = Matrix("A", (n, n))
        A.set_property(properties.FULL_RANK)

        B = Matrix("B", (n, n))
        B.set_property(properties.FULL_RANK)

        C = Matrix("C", (n, n))
        C.set_property(properties.FULL_RANK)

        D = Matrix("D", (n, n))
        D.set_property(properties.FULL_RANK)

        alpha = Scalar("alpha")

        X = Matrix("X", (n, n))

        self.eqns = Equations(
                        Equal(X, Plus(Transpose(A), B, Times(alpha, C), Transpose(D)))
                        )


class Example161():
    def __init__(self):

        # TAGS
        # SAC

        n = 5000

        A = Matrix("A", (n, n))
        A.set_property(properties.FULL_RANK)

        B = Matrix("B", (n, n))
        B.set_property(properties.FULL_RANK)

        # C = Matrix("C", (n, n))
        # C.set_property(properties.FULL_RANK)

        D = Matrix("D", (n, n))
        D.set_property(properties.FULL_RANK)
        D.set_property(properties.DIAGONAL)

        alpha = Scalar("alpha")

        X = Matrix("X", (n, n))

        self.eqns = Equations(
                        Equal(X, Plus(A, Times(D, B)))
                        )


class Example162():
    def __init__(self):

        # TAGS
        # linear system, simple, ldlt

        n = 2000
        m = 10

        A = Matrix("A", (n, m))

        D = Matrix("D", (n, n))
        D.set_property(properties.SYMMETRIC)
        D.set_property(properties.FULL_RANK)

        W = Matrix("W", (n, m))

        self.eqns = Equations(Equal(W, Times(Inverse(D), A)))


class Example163():
    def __init__(self):

        # TAGS
        # 

        n = 1000

        A = Matrix("A", (n, n))
        A.set_property(properties.FULL_RANK)

        B = Matrix("B", (n, n))
        B.set_property(properties.FULL_RANK)

        C = Matrix("C", (n, n))
        C.set_property(properties.FULL_RANK)

        D = Matrix("D", (n, n))
        D.set_property(properties.FULL_RANK)


        X = Matrix("X", (n, n))
        Y = Matrix("Y", (n, n))
        Z = Matrix("Z", (n, n))

        self.eqns = Equations(
                        Equal(X, Transpose(A)),
                        Equal(Y, Times(B, Plus(C, D)))
                        )