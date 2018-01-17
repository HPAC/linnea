from linnea.algebra.expression import Matrix, Vector, Equal, Times, Inverse, Transpose
from linnea.algebra.equations import Equations
from linnea.algebra.properties import Property

n = 1500
m = 1000

X = Matrix("X", (n, m))
X.set_property(Property.INPUT)
X.set_property(Property.FULL_RANK)

y = Vector("y", (n, 1))
y.set_property(Property.INPUT)

b = Vector("b", (m, 1))
b.set_property(Property.OUTPUT)

# b = (X^T X)^-1 X^T y
equations = Equations(Equal(b, Times(Inverse(Times(Transpose(X), X)), Transpose(X), y)))