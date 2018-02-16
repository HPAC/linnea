from ..algebra.expression import Scalar, Matrix, \
                                Equal, Plus, Times, Transpose, Inverse, \
                                Symbol, ConstantScalar

from ..algebra.properties import Property as properties
from ..algebra.transformations import simplify
from ..algebra.representations import to_SOP

from .graph.base import base

from .. import config
from .. import temporaries

from .matrix_sum import decompose_sum

import matchpy
import copy

collections_module = config.import_collections()

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
                    SYM1.has_property(properties.ORTHOGONAL) and
                    SYM2.has_property(properties.DIAGONAL) and
                    WD1.has_property(properties.SCALAR) and
                    SYM3.has_property(properties.IDENTITY)
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
    tmp = temporaries.create_tmp(to_SOP(simplify(diagonal_sum)), True)
    
    if substitution["WS1"]:
        replacement = Plus(Times(Transpose(substitution["SYM1"]), tmp, substitution["SYM1"]), *substitution["WS1"])
    else:
        replacement = Times(Transpose(substitution["SYM1"]), tmp, substitution["SYM1"])
    equations_list[eqn_idx] = matchpy.replace(equations_list[eqn_idx], (1,)+tuple(position), replacement)
    equations_list[eqn_idx] = simplify(equations_list[eqn_idx])

    new_equation = Equal(tmp, diagonal_sum)
    equations_list.insert(eqn_idx, new_equation)
    new_equations = Equations(*equations_list)

    return (new_equations, new_equations.metric(), base.EdgeLabel())

eigen2 = matchpy.Pattern(
            Plus(Times(SYM1, SYM2, Transpose(SYM1)), Times(WD1, SYM3), WS1),
            matchpy.CustomConstraint(
                lambda SYM1, SYM2, SYM3, WD1:
                    SYM1.has_property(properties.ORTHOGONAL) and
                    SYM2.has_property(properties.DIAGONAL) and
                    WD1.has_property(properties.SCALAR) and
                    SYM3.has_property(properties.IDENTITY)
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
    tmp = temporaries.create_tmp(to_SOP(simplify(diagonal_sum)), True)

    if substitution["WS1"]:
        replacement = Plus(Times(substitution["SYM1"], tmp, Transpose(substitution["SYM1"])), *substitution["WS1"])
    else:
        replacement = Times(substitution["SYM1"], tmp, Transpose(substitution["SYM1"]))
    equations_list[eqn_idx] = matchpy.replace(equations_list[eqn_idx], (1,)+tuple(position), replacement)
    equations_list[eqn_idx] = simplify(equations_list[eqn_idx])

    new_equation = Equal(tmp, diagonal_sum)
    equations_list.insert(eqn_idx, new_equation)
    new_equations = Equations(*equations_list)

    return (new_equations, new_equations.metric(), base.EdgeLabel())

def symmetric_product_constraint(WP1, WP2, _A):
    p1 = Times(*WP1)
    p2 = Times(*WP2)
    return _A.has_property(properties.ADMITS_FACTORIZATION) and _A.has_property(properties.SPD) and p1.transpose_of(p2)

symmetric_product = matchpy.Pattern(
                        Times(WP1, _A, WP2),
                        matchpy.CustomConstraint(symmetric_product_constraint)
                    )

def symmetric_product_callback(substitution, equations, eqn_idx, position):
    # symmetric product
    equations_list = list(equations.equations)

    matched_kernel = collections_module.cholesky.set_match({"_A": substitution["_A"]}, False, CSE_rules=True)

    replacement = Times(*substitution["WP1"], matched_kernel.replacement, *substitution["WP2"])

    equations_list[eqn_idx] = matchpy.replace(equations_list[eqn_idx], (1,)+tuple(position), replacement)

    # deal with additional occurrences of the replaced subexpression
    # common_subexp_rules = cholesky.get_rules(substitution)
    # common_subexp_rules = cholesky.get_CSE_rules()
    # common_subexp_rules = matched_kernel.CSE_rules

    new_equations = Equations(*equations_list)
    new_equations = new_equations.replace_all(matched_kernel.CSE_rules)
    new_equations.set_equivalent(equations)
    
    edge_label = base.EdgeLabel(matched_kernel)
    return (new_equations, new_equations.metric(), edge_label)


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

    equations_list[eqn_idx] = matchpy.replace(equations_list[eqn_idx], [1]+list(position), replacement)
    
    new_equations = Equations(*equations_list)

    edge_label = base.EdgeLabel(*matched_kernels)
    return (new_equations, new_equations.metric(), edge_label)


# A^T B + B^T A + A^T C A (C is symmetric)
trick4 = matchpy.Pattern(
            Plus(Times(Transpose(WD1), WD2), Times(Transpose(WD2), WD1), Times(Transpose(WD1), WD3, WD1), WS1),
            matchpy.CustomConstraint(
                lambda WD1, WD2, WD3, WS1: WD3.has_property(properties.SYMMETRIC)
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
    tmp = temporaries.create_tmp(to_SOP(simplify(sum_expr)), True)
    new_equation = Equal(tmp, sum_expr)

    replacement = Plus(Times(Transpose(substitution["WD1"]), tmp), Times(Transpose(tmp), substitution["WD1"]), *substitution["WS1"])

    equations_list[eqn_idx] = matchpy.replace(equations_list[eqn_idx], [1]+list(position), replacement)

    equations_list = equations_list.insert(eqn_idx, new_equation)
    new_equations = Equations(*equations_list)

    return (new_equations, new_equations.metric(), base.EdgeLabel())

# patterns = [eigen1, eigen2, symmetric_product, trick3]

trick_MA = matchpy.ManyToOneMatcher()
trick_MA.add(eigen1, label=eigen1_callback)
trick_MA.add(eigen2, label=eigen2_callback)
trick_MA.add(symmetric_product, label=symmetric_product_callback)
trick_MA.add(trick3, label=trick3_callback)
trick_MA.add(trick4, label=trick4_callback)

# automaton = PMA()
# automaton.add_patterns(patterns)
# automaton.to_dot_file()