from linnea.algebra.expression import Matrix, Vector, Equal, Times, Inverse
from linnea.algebra.equations import Equations
from linnea.algebra.properties import Property

n = 2000

A = Matrix("A", (n, n))
B = Matrix("B", (n, n))
B.set_property(Property.SPD)
c = Vector("c", (n, 1))
x = Matrix("x", (n, 1))

equations = Equations(
                    Equal(x, Times(A, Inverse(B), c))
                    )