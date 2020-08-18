from linnea.algebra.expression import Matrix, Vector, Scalar, Plus, Times, \
                                      Inverse, Transpose, InverseTranspose, \
                                      IdentityMatrix, ZeroMatrix, ConstantScalar
from linnea.algebra.transformations import simplify
from linnea.algebra.properties import Property

import pytest

M1 = Matrix("M1", (100, 100))
M2 = Matrix("M2", (100, 100))
M3 = Matrix("M3", (100, 100))
M4 = Matrix("M3", (200, 200))
S = Matrix("S", (100, 100), properties={Property.SYMMETRIC})
SPD = Matrix("SPD", (100, 100), properties={Property.SPD})
Q1 = Matrix("Q1", (100, 100), properties={Property.ORTHOGONAL})
Q2 = Matrix("Q2", (100, 100), properties={Property.ORTHOGONAL})
I = IdentityMatrix(100, 100)
Ic = IdentityMatrix(200, 100)
Ir = IdentityMatrix(100, 200)
Zero = ZeroMatrix(100, 100)
Qc = Matrix("Qc", (200, 100), properties={Property.ORTHOGONAL_COLUMNS})
Qr = Matrix("Qr", (100, 200), properties={Property.ORTHOGONAL_ROWS})
alpha = Scalar("alpha")
beta = Scalar("beta")
beta = Scalar("gamma")
zero = ConstantScalar(0.0)
one = ConstantScalar(1.0)
two = ConstantScalar(2.0)
four = ConstantScalar(4.0)
pointfive = ConstantScalar(0.5)
minusone = ConstantScalar(-1.0)

@pytest.mark.parametrize(
    "input,                                 expected_output",
    [
        (Transpose(S),                      S),
        # inverse
        (Times(M1, Inverse(M1)),            I),
        (Times(Inverse(M1), M1),            I),
        (Times(M2, M1, Inverse(M1)),        M2),
        (Times(M2, Inverse(M1), M1),        M2),
        (Times(M1, Inverse(M1), M2),        M2),
        (Times(Inverse(M1), M1, M2),        M2),
        (Times(M2, M1, Inverse(M1), M3),    Times(M2, M3)),
        (Times(M2, Inverse(M1), M1, M3),    Times(M2, M3)),
        # TODO Should this be simplified? What if alpha is zero?
        # (Times(alpha, Inverse(alpha)),          one),
        # (Times(Inverse(alpha), alpha),          one),
        # (Times(alpha, beta, Inverse(alpha)),    beta),
        # (Times(Inverse(alpha), beta, alpha),    beta),
        # (Times(alpha, M1, Inverse(alpha)),      M1),
        # (Times(Inverse(alpha), M1, alpha),      M1),
        # orthogonal
        (Inverse(Q1),                       Transpose(Q1)),
        (Times(Q1, Transpose(Q1)),          I),
        (Times(Transpose(Q1), Q1),          I),
        (Times(M1, Q1, Transpose(Q1)),      M1),
        (Times(M1, Transpose(Q1), Q1),      M1),
        (Times(Q1, Transpose(Q1), M1),      M1),
        (Times(Transpose(Q1), Q1, M1),      M1),
        (Times(M1, Q1, Transpose(Q1), M2),  Times(M1, M2)),
        (Times(M1, Transpose(Q1), Q1,M2),   Times(M1, M2)),
        (Times(Q1, Transpose(Q2)),          Times(Q1, Transpose(Q2))),
        (Times(Transpose(Q2), Q1),          Times(Transpose(Q2), Q1)),
        # orthogonal columns
        (Times(Transpose(Qc), Qc),          I),
        (Times(M1, Transpose(Qc), Qc),      M1),
        (Times(Transpose(Qc), Qc, M1),      M1),
        (Times(M1, Transpose(Qc), Qc,M2),   Times(M1, M2)),
        (Times(Qc, Transpose(Qc)),          Times(Qc, Transpose(Qc))),
        # orthogonal rows
        (Times(Qr, Transpose(Qr)),          I),
        (Times(M1, Qr, Transpose(Qr)),      M1),
        (Times(Qr, Transpose(Qr), M1),      M1),
        (Times(M1, Qr, Transpose(Qr), M2),  Times(M1, M2)),
        (Times(Transpose(Qr), Qr),          Times(Transpose(Qr), Qr)),
        # Identity
        (Transpose(I),                      I),
        (Inverse(I),                        I),
        (InverseTranspose(I),               I),
        (Times(I, I),                       I),
        (Times(M1, I),                      M1),
        (Times(I, M1),                      M1),
        (Times(M1, I, M2),                  Times(M1, M2)),
        (Times(M1, Ir, M4),                 Times(M1, Ir, M4)),
        (Times(M4, Ic, M1),                 Times(M4, Ic, M1)),
        (Times(M1, I, I, M2),               Times(M1, M2)),
        # Zero
        (Transpose(Zero),                   Zero),
        (Plus(Zero, Zero),                  Zero),
        (Plus(M1, Zero),                    M1),
        (Plus(Zero, M1),                    M1),
        (Plus(M1, Zero, M2),                Plus(M1, M2)),
        (Plus(M1, Zero, Zero, M2),          Plus(M1, M2)),
        (Times(M1, Zero),                   Zero),
        (Times(Zero, M1),                   Zero),
        (Times(M1, Zero, M2),               Zero),
        # scalars
        (Times(two, two),                   four),
        (Plus(two, two),                    four),
        (Plus(zero, one),                   one),
        (Times(two, pointfive),             one),
        (Times(one, M1),                    M1),
        (Times(two, pointfive, M1),         M1),
        (Times(two, M1, pointfive),         M1),
        (Times(M1, two, pointfive),         M1),
        (Times(zero, M1),                   Zero),
        (Inverse(two),                      pointfive),
        (Transpose(two),                    two),
        (InverseTranspose(two),             pointfive),
        (Times(Inverse(two), two),          one),
        (Times(two, Inverse(two)),          one),
        (Plus(two, minusone),               one),
        (Plus(one, minusone),               zero),
        # addition
        (Plus(M1, M1),                          Times(two, M1)),
        (Plus(M1, Times(minusone, M1)),         Zero),
        (Plus(alpha, Times(minusone, alpha)),   zero),
        (Plus(M1, M1, M2),                      Plus(Times(two, M1), M2)),
        (Plus(M1, Times(minusone, M1), M2),     M2),
        # sorting
        (Times(two, alpha),                  Times(two, alpha)),
        (Times(alpha, two),                  Times(two, alpha)),
        (Times(alpha, M1),                   Times(alpha, M1)),
        (Times(M1, alpha),                   Times(alpha, M1)),
        (Times(two, M1),                     Times(two, M1)),
        (Times(M1, two),                     Times(two, M1)),
        (Times(two, alpha, M1),              Times(two, alpha, M1)),
        (Times(alpha, two, M1),              Times(two, alpha, M1)),
        (Times(alpha, M1, two),              Times(two, alpha, M1)),
        (Times(two, M1, alpha),              Times(two, alpha, M1)),
        (Times(M1, two, alpha),              Times(two, alpha, M1)),
        (Times(M1, alpha, two),              Times(two, alpha, M1)),
    ]
)
def test_simplify(input, expected_output):
    output = simplify(input)
    assert output == expected_output, "{} did not simplify to {}, but instead to {}.".format(input, expected_output, output)

