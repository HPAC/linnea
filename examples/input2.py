from linnea.algebra.expression import Matrix, Equal, Times, Inverse, Transpose, InverseTranspose
from linnea.algebra.equations import Equations
from linnea.algebra.properties import Property

n = 100

A = Matrix("A", (n, n))
A.set_property(Property.INPUT)

B = Matrix("B", (n, n))
B.set_property(Property.INPUT)
B.set_property(Property.SPD)

C = Matrix("C", (n, n))
C.set_property(Property.INPUT)

D = Matrix("D", (n, n))
D.set_property(Property.INPUT)

X = Matrix("X", (n, n))
X.set_property(Property.OUTPUT)

Y = Matrix("Y", (n, n))
Y.set_property(Property.OUTPUT)

equations = Equations(
                    Equal(X, Times(A, InverseTranspose(B), C)),
                    Equal(Y, Times(Inverse(B), Transpose(A), D)),
                    )