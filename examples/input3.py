from linnea.algebra.expression import Matrix, Vector, Equal, Times, Inverse
from linnea.algebra.equations import Equations
from linnea.algebra.properties import Property

n = 2000

A = Matrix("A", (n, n))
A.set_property(Property.INPUT)

B = Matrix("B", (n, n))
B.set_property(Property.INPUT)
B.set_property(Property.SPD)

c = Vector("c", (n, 1))
c.set_property(Property.INPUT)

x = Matrix("x", (n, 1))
x.set_property(Property.OUTPUT)

equations = Equations(
                    Equal(x, Times(A, Inverse(B), c))
                    )