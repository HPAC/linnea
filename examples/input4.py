from linnea.algebra.expression import Matrix, Vector, Equal, Times, Inverse, Transpose
from linnea.algebra.equations import Equations
from linnea.algebra.properties import Property

n = 1500
m = 1000

X = Matrix("X", (n, m))
X.set_property(Property.FULL_RANK)
M = Matrix("M", (n, n))
M.set_property(Property.SPD)
y = Vector("y", (n, 1))
b = Vector("b", (m, 1))

equations = Equations(Equal(b, Times(Inverse(Times(Transpose(X), Inverse(M), X ) ), Transpose(X), Inverse(M), y)))
