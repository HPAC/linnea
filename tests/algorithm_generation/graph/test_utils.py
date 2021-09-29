from linnea.algebra.expression import Matrix, Vector, Scalar, Plus, Times, \
                                      Inverse, Transpose, InverseTranspose, \
                                      IdentityMatrix, ZeroMatrix, ConstantScalar
from linnea.algorithm_generation.graph.utils import is_simple_summand, is_simple_plus

import pytest

M1 = Matrix("M1", (100, 100))
M2 = Matrix("M2", (100, 100))
M3 = Matrix("M3", (100, 100))
v1 = Matrix("v1", (100, 100))
v2 = Matrix("v2", (100, 100))
alpha = Scalar("alpha")
two = ConstantScalar(2.0)


@pytest.mark.parametrize(
        "input,                                 expected_output",
    [
        (M1,                                    False),
        (v1,                                    False),
        (alpha,                                 False),
        (two,                                   False),
        (Plus(M1, M2),                          True),
        (Plus(M1, Times(alpha, M2)),            True),
        (Plus(M1, Times(two, M2)),              True),
        (Plus(M1, Times(M2, M3)),               False),
        (Plus(M1, M2, M3),                      True),
        (Plus(M1, Times(alpha, M2), M3),        True),
        (Plus(M1, Times(two, M2), M3),          True),
    ]
)
def test_is_simple_plus(input, expected_output):
    output = is_simple_plus(input)
    assert output == expected_output, "is_simple_plus returned {} for {}, but should have returned {}.".format(output, input, expected_output)


@pytest.mark.parametrize(
        "input,                                 expected_output",
    [
        (M1,                                    True),
        (Transpose(M1),                         True),
        (Inverse(M1),                           False),
        (InverseTranspose(M1),                  False),
        (v1,                                    True),
        (Transpose(v1),                         True),
        (alpha,                                 True),
        (two,                                   True),
        (Times(alpha, M1),                      True),
        (Times(alpha, Transpose(M1)),           True),
        (Times(alpha, Inverse(M1)),             False),
        (Times(alpha, InverseTranspose(M1)),    False),
        (Times(alpha, v1),                      True),
        (Times(alpha, Transpose(v1)),           True),
        (Times(two, M1),                        True),
        (Times(two, Transpose(M1)),             True),
        (Times(two, Inverse(M1)),               False),
        (Times(two, InverseTranspose(M1)),      False),
        (Times(two, alpha),                     True),
        (Times(M1, M2),                         False),
        (Times(M1, v1),                         False),
        (Times(two, alpha, M1),                 False),
        (Times(two, M1, M2),                    False),
        (Times(alpha, M1, M2),                  False),
        (Times(M1, Plus(M2, M3)),               False),
        (Times(alpha, Plus(M2, M3)),            False),
        (Times(two, Plus(M2, M3)),              False),
        (Times(alpha, Plus(v1, v2)),            False),
        (Times(two, Plus(v1, v2)),              False),
    ]
)
def test_is_simple_summand(input, expected_output):
    output = is_simple_summand(input)
    assert output == expected_output, "is_simple_summand returned {} for {}, but should have returned {}.".format(output, input, expected_output)


if __name__ == '__main__':
    pytest.main()