from ..algebra.expression import Scalar, Matrix, \
                                Equal, Plus, Times, Transpose, Inverse, \
                                Symbol, ConstantScalar

from ..algebra.properties import Property
from ..algebra.transformations import simplify
from ..algebra.representations import to_SOP
from ..algebra.equations import Equations

from .. import config
from .. import temporaries

from ..utils import is_inverse

from .reductions import decompose_sum

import matchpy
import copy

collections_module = config.import_collections()

def apply_tricks(equations):
    equation, eqn_idx = equations.process_next()
    if not equation:
        return
    for expr, position in equation.rhs.preorder_iter():
        for callback_func, substitution in trick_MA.match(expr):
            result = callback_func(substitution, equations, eqn_idx, position)
            if result:
                yield result


WD1 = matchpy.Wildcard.dot("WD1")
WD2 = matchpy.Wildcard.dot("WD2")
WD3 = matchpy.Wildcard.dot("WD3")
SYM1 = matchpy.Wildcard.symbol("SYM1")
SYM2 = matchpy.Wildcard.symbol("SYM2")
SYM3 = matchpy.Wildcard.symbol("SYM3")
_A = matchpy.Wildcard.symbol("_A") # It's important that this Wildcard has the same name as the one in the pattern for Cholesky
WS1 = matchpy.Wildcard.star("WS1")
WP1 = matchpy.Wildcard.plus("WP1")
WP2 = matchpy.Wildcard.plus("WP2")

eigen1 = matchpy.Pattern(
            Plus(Times(Transpose(SYM1), SYM2, SYM1), Times(WD1, SYM3), WS1),
            matchpy.CustomConstraint(
                lambda SYM1, SYM2, SYM3, WD1:
                    SYM1.has_property(Property.ORTHOGONAL) and
                    SYM2.has_property(Property.DIAGONAL) and
                    WD1.has_property(Property.SCALAR) and
                    SYM3.has_property(Property.IDENTITY)
            )
        )

def eigen1_callback(substitution, equations, eqn_idx, position):
    # Here, the "Eigen-Trick" is applied. That is, given
    # Plus([Q^T W Q + alpha I]), tmp = W + alpha I is extraced and the
    # entire expression is replaced with Times([Q^T tmp Q])
    # The sum can not be computed directly with decompose_sum
    # because it is not necessarily a sufficiently simple
    # sum.
    equations_list = list(equations.equations)

    diagonal_sum = Plus(substitution["SYM2"], Times(substitution["WD1"], substitution["SYM3"]))
    tmp = temporaries.create_tmp(diagonal_sum, True)
    
    if substitution["WS1"]:
        replacement = Plus(Times(Transpose(substitution["SYM1"]), tmp, substitution["SYM1"]), *substitution["WS1"])
    else:
        replacement = Times(Transpose(substitution["SYM1"]), tmp, substitution["SYM1"])
    equations_list[eqn_idx] = matchpy.replace(equations_list[eqn_idx], (1,)+tuple(position), replacement)
    equations_list[eqn_idx] = simplify(equations_list[eqn_idx])

    new_equation = Equal(tmp, diagonal_sum)
    equations_list.insert(eqn_idx, new_equation)
    new_equations = Equations(*equations_list)

    return (new_equations, ())

eigen2 = matchpy.Pattern(
            Plus(Times(SYM1, SYM2, Transpose(SYM1)), Times(WD1, SYM3), WS1),
            matchpy.CustomConstraint(
                lambda SYM1, SYM2, SYM3, WD1:
                    SYM1.has_property(Property.ORTHOGONAL) and
                    SYM2.has_property(Property.DIAGONAL) and
                    WD1.has_property(Property.SCALAR) and
                    SYM3.has_property(Property.IDENTITY)
            )
        )


def eigen2_callback(substitution, equations, eqn_idx, position):
    # Here, the "Eigen-Trick" is applied. That is, given
    # Plus([Q W Q^T + alpha I]), tmp = W + alpha I is extraced and the
    # entire expression is replaced with Times([Q tmp Q^T])
    # The sum can not be computed directly with decompose_sum
    # because it is not necessarily a sufficiently simple
    # sum.
    equations_list = list(equations.equations)

    diagonal_sum = Plus(substitution["SYM2"], Times(substitution["WD1"], substitution["SYM3"]))
    tmp = temporaries.create_tmp(diagonal_sum, True)

    if substitution["WS1"]:
        replacement = Plus(Times(substitution["SYM1"], tmp, Transpose(substitution["SYM1"])), *substitution["WS1"])
    else:
        replacement = Times(substitution["SYM1"], tmp, Transpose(substitution["SYM1"]))
    equations_list[eqn_idx] = matchpy.replace(equations_list[eqn_idx], (1,)+tuple(position), replacement)
    equations_list[eqn_idx] = simplify(equations_list[eqn_idx])

    new_equation = Equal(tmp, diagonal_sum)
    equations_list.insert(eqn_idx, new_equation)
    new_equations = Equations(*equations_list)

    return (new_equations, ())

def symmetric_product_constraint(WP1, WP2, _A):
    p1 = Times(*WP1)
    p2 = Times(*WP2)
    return _A.has_property(Property.ADMITS_FACTORIZATION) and _A.has_property(Property.SPD) and p1.transpose_of(p2)

symmetric_product = matchpy.Pattern(
                        Times(WP1, _A, WP2),
                        matchpy.CustomConstraint(symmetric_product_constraint)
                    )

def symmetric_product_callback(substitution, equations, eqn_idx, position):
    # symmetric product
    equations_list = list(equations.equations)

    # This trick is not applied if the current position is inside an inverse
    # because Cholesky is applied in that case anyway.
    expr = equations_list[eqn_idx].rhs
    for p in position:
        if is_inverse(expr):
            return
        expr = expr[p]
    # There is no need to check the type of expr here again because it is Times.

    matched_kernel = collections_module.cholesky.set_match({"_A": substitution["_A"]}, False)

    replacement = Times(*substitution["WP1"], matched_kernel.replacement, *substitution["WP2"])

    equations_list[eqn_idx] = matchpy.replace(equations_list[eqn_idx], (1,)+tuple(position), replacement)

    new_equations = Equations(*equations_list)
    new_equations.set_equivalent(equations)
    
    return (new_equations, (matched_kernel,))


# A^T B + B^T A + A^T A
trick3 = matchpy.Pattern(Plus(Times(Transpose(WD1), WD2), Times(Transpose(WD2), WD1), Times(Transpose(WD1), WD1), WS1))


def trick3_callback(substitution, equations, eqn_idx, position):
    # A^T B + B^T A + A^T A
    # A^T (B + 1/2 A) + (B + 1/2 A)^T A
    # WD1 = A
    # WD2 = B

    equations_list = list(equations.equations)
    one_half = ConstantScalar(0.5)
    sum_expr = Plus(Times(one_half, substitution["WD1"]), substitution["WD2"])
    tmp, matched_kernels = decompose_sum(sum_expr)

    replacement = Plus(Times(Transpose(substitution["WD1"]), tmp), Times(Transpose(tmp), substitution["WD1"]), *substitution["WS1"])

    equations_list[eqn_idx] = matchpy.replace(equations_list[eqn_idx], (1,)+position, replacement)
    
    new_equations = Equations(*equations_list)

    return (new_equations, matched_kernels)


# A^T B + B^T A + A^T C A (C is symmetric)
trick4 = matchpy.Pattern(
            Plus(Times(Transpose(WD1), WD2), Times(Transpose(WD2), WD1), Times(Transpose(WD1), WD3, WD1), WS1),
            matchpy.CustomConstraint(
                lambda WD1, WD2, WD3, WS1: WD3.has_property(Property.SYMMETRIC)
            )
        )

def trick4_callback(substitution, equations, eqn_idx, position):
    # A^T B + B^T A + A^T C A (C is symmetric)
    # A^T (B + 1/2 C A) + (B + 1/2 C A)^T A
    # WD1 = A
    # WD2 = B
    # WD3 = C

    equations_list = list(equations.equations)

    one_half = ConstantScalar(0.5)
    sum_expr = Plus(Times(one_half, substitution["WD3"], substitution["WD1"]), substitution["WD2"])
    tmp = temporaries.create_tmp(sum_expr, True)
    new_equation = Equal(tmp, sum_expr)

    replacement = Plus(Times(Transpose(substitution["WD1"]), tmp), Times(Transpose(tmp), substitution["WD1"]), *substitution["WS1"])

    equations_list[eqn_idx] = matchpy.replace(equations_list[eqn_idx], (1,)+position, replacement)

    equations_list.insert(eqn_idx, new_equation)
    new_equations = Equations(*equations_list)

    return (new_equations, ())

def syr2k1_constraint(WP1, WP2, _A):
    p1 = Times(*WP1)
    p2 = Times(*WP2)
    return _A.has_property(Property.SYMMETRIC) and p1.transpose_of(p2)

syr2k1 = matchpy.Pattern(
                        Times(WP1, _A, WP2),
                        matchpy.CustomConstraint(syr2k1_constraint)
                    )

def syr2k1_callback(substitution, equations, eqn_idx, position):
    # A B A^T (B is symmetric)
    # 1/2 A (A B)^T + 1/2 (A B) A^T
    # TODO what about A^T B A? Is it necessary to treat it differently? 
    equations_list = list(equations.equations)

    one_half = ConstantScalar(0.5)
    new_expr = Times(*substitution["WP1"], substitution["_A"])
    tmp = temporaries.create_tmp(new_expr, True)
    new_equation = Equal(tmp, new_expr)

    A = Times(*substitution["WP1"])
    replacement = Plus(Times(one_half, A, Transpose(tmp)), Times(one_half, tmp, Transpose(A)))

    equations_list[eqn_idx] = matchpy.replace(equations_list[eqn_idx], (1,)+position, replacement)

    equations_list.insert(eqn_idx, new_equation)
    new_equations = Equations(*equations_list)

    return (new_equations, ())

# patterns = [eigen1, eigen2, symmetric_product, trick3]

trick_MA = matchpy.ManyToOneMatcher()
trick_MA.add(eigen1, label=eigen1_callback)
trick_MA.add(eigen2, label=eigen2_callback)
trick_MA.add(symmetric_product, label=symmetric_product_callback)
trick_MA.add(trick3, label=trick3_callback)
trick_MA.add(trick4, label=trick4_callback)
trick_MA.add(syr2k1, label=syr2k1_callback)

# automaton = PMA()
# automaton.add_patterns(patterns)
# automaton.to_dot_file()