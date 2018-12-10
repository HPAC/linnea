from linnea.algebra.expression import Matrix, Equal, Times, Inverse, Transpose, InverseTranspose
from linnea.algebra.equations import Equations
from linnea.algebra.properties import Property

n = 100

A = Matrix("A", (n, n))
B = Matrix("B", (n, n))
B.set_property(Property.SPD)
C = Matrix("C", (n, n))
D = Matrix("D", (n, n))
X = Matrix("X", (n, n))
Y = Matrix("Y", (n, n))


equations = Equations(
                    Equal(X, Times(A, InverseTranspose(B), C)),
                    Equal(Y, Times(Inverse(B), Transpose(A), D)),
                    )