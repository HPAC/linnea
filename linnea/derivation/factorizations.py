from ..algebra.expression import Symbol, Times, Plus
from ..algebra.equations import Equations
from ..utils import powerset, is_inverse, is_transpose

from .. import config

from enum import Enum, unique
from collections import namedtuple

import itertools
import matchpy

collections_module = config.import_collections()

@unique
class InverseType(Enum):
    linear_system = 1
    explicit_inversion = 2
    none = 3

Occurrence = namedtuple('Occurrence', ['eqn_idx', 'position', 'operand', 'type', 'group', 'symbol'])

def apply_factorizations(equations, operands_to_factor, factorization_dict):
    """This function generates new equations by applying factorizations.

    For applying factorizations, a number of rules are applied:
    - Matrices are only factored if they have the ADMITS_FACTORIZATION
      property (this is tested by find_operands_to_factor() ).
    - We never apply different factorizations to different occurrences of
      the same matrix. As an example, if the matrix A shows up twice, we
      will never apply LU to one occurrence and QR to another.
    - Factorization will always be applied to matrices that appear
      immediately inside an inverse. That is, the A in Inverse(A) will be
      factored. A and B in Inverse(Times(A, B)) don't have to be factored.
    - If there is a summand that contains a matrix which is factored, but
      the summand does not contain any occurrences of that matrix within an
      inverse, that matrix will not be factored in this summand. This is
      to make sure that for A+inv(A), the first A is not factored.
    - Some factorizations rule out others. If Cholesky can be applied, LU
      will no be applied. Which factorization can be applied per operand is
      decided in DS_factorizations(). The factorization_dict contains the
      valid factorizations.

    This function applies factorization in all possible combinations that
    obey the rules above.
    """

    # find all occurrences
    all_occurrences = list(find_occurrences(equations, operands_to_factor))

    blocking_products = list(find_blocking_products(equations, operands_to_factor))

    # Removing groups (summands) which do not contain any inverted occurrences.
    candidate_occurrences = []
    for oc_group in group_occurrences(all_occurrences):
        if any(oc.type != InverseType.none for oc in oc_group):
            candidate_occurrences.extend(oc_group)

    # collect all operands that show up
    ops = set(oc.operand.name for oc in candidate_occurrences)

    # Symbols directely inside an inverse always have to be factored.
    ops_must_factor = set()
    ops_may_factor = set()
    for op in ops:
        if any(oc.operand.name == op and oc.symbol for oc in candidate_occurrences):
            ops_must_factor.add(op)
        else:
            ops_may_factor.add(op)

    # sorting here removes randomness
    for ops_subset in powerset(sorted(ops_may_factor)):

        factor_ops = ops_must_factor.union(ops_subset)
        if not factor_ops or factor_ops in blocking_products:
            continue

        factorizations_candidates = []
        # sorting here removes randomness
        factor_ops_sorted = sorted(factor_ops)
        for op in factor_ops_sorted:
            factorizations_candidates.append(factorization_dict[op])

        # apply all factorizations
        for factorizations in itertools.product(*factorizations_candidates):
            facts_dict = dict(zip(factor_ops_sorted, factorizations))

            # collect matched kernels (avoiding duplicates)
            matched_kernels = []
            _already_seen = set()
            for matched_kernel in factorizations:
                if matched_kernel.id not in _already_seen:
                    matched_kernels.append(matched_kernel)
                    _already_seen.add(matched_kernel.id)

            # collect replacements 
            replacements_per_equation = dict()

            for oc in candidate_occurrences:
                if oc.operand.name in factor_ops:
                    replacements_per_equation.setdefault(oc.eqn_idx, []).append((oc.position, facts_dict[oc.operand.name].replacement))

            # replace
            equations_list = list(equations.equations)

            for eqn_idx, replacements in replacements_per_equation.items():
                if replacements:
                    equations_list[eqn_idx] = matchpy.replace_many(equations_list[eqn_idx], replacements)
            
            equations_copy = Equations(*equations_list)
            equations_copy = equations_copy.simplify()
            equations_copy.set_equivalent(equations)
            equations_copy = equations_copy.to_SOP().simplify()

            yield (equations_copy, matched_kernels)


def construct_factorization_dict(operands_to_factor):
    
    factorization_dict = dict()
    for operand in operands_to_factor:
        for type in collections_module.factorizations_by_type:
            for kernel in type:
                matches = list(matchpy.match(operand, kernel.pattern))

                # This is important. Otherwise, the "break" at the end of
                # the loop is reached. In that case, only the first
                # factorization of each type is applied.
                if not matches:
                    continue

                # this is kind of stupid, there is only one match
                for match_dict in matches:
                    matched_kernel = kernel.set_match(match_dict, False)
                    factorization_dict.setdefault(operand.name, []).append(matched_kernel)
                break
    return factorization_dict


def find_occurrences(equations, operands_to_factor):
    """Finds all occurrences of operands that have to be factored.

    Finds all occurrences of operands that have to be factored. This includes
    occurrences which are not in inverses.

    This function returns and Occurrences object, which has the following
    attributes:
    * eqn_idx: The index of the equation of the occurrence.
    * position: The position of the occurrence in equation.
    * operand: The Operand.
    * type: The type of this occurrence. This is an InverseType object.
    * group: An identifier for the group. Either the root of the expression, or
      the position of this subexpression in the last sum.
    * symbol: True if this occurrence is a symbol inverse, that is Inverse(A).

    Args:
        equations (Equations): The equations that are searched.
        operands_to_factor (set): Set of operands to search for.

    Yields:
        Occurrence: All occurrences of the operands.

    """

    for eqn_idx, equation in enumerate(equations):
        for res in _find_occurrences(equation.rhs, operands_to_factor, inv_type=InverseType.none, position=(1,), group=(1,)):
            # for grouping, we also need the eqn_idx
            yield Occurrence(eqn_idx, *res)


def _find_occurrences(expr, operands_to_factor, inv_type=InverseType.none, position=(), group=None, symbol=False, predecessor=None):

    if isinstance(expr, Symbol):
        if expr in operands_to_factor:
            yield (position, expr, inv_type, group, symbol)
        return

    if is_inverse(expr):
        if isinstance(predecessor, Times):
            inv_type = InverseType.linear_system
        else:
            inv_type = InverseType.explicit_inversion
        symbol = isinstance(expr.operand, Symbol)

    for n, operand in enumerate(expr.operands):
        new_position = position + (n,)
        new_group = group
        if isinstance(expr, Plus):
            new_group = new_position
        yield from _find_occurrences(operand, operands_to_factor, inv_type, new_position, new_group, symbol, expr)


def group_occurrences(occurrences):
    # group symbol inverses by eqn_idx, inverse group, operand
    occurrences = sorted(occurrences, key=grouping_keyfunc)
    occurrences_grouped = []
    # grouping by eqn_idx, inverse group, operand
    for key, group in itertools.groupby(occurrences, grouping_keyfunc):
        occurrences_grouped.append(list(group))
    return occurrences_grouped


def grouping_keyfunc(oc):
    # eqn_idx, operand, group
    if oc.group is None:
        return (oc.eqn_idx, oc.operand, [])
    else:
        return (oc.eqn_idx, oc.operand, oc.group)


def find_blocking_products(equations, operands_to_factor):
    """Identifies sets of operands likely to lead to dead ends when factored.

    In many cases, when all operands in a product are factored, the resulting
    expression can not be computed anymore and becomes a dead end. The purpose
    of this function is to identify such products.
    """

    for equation in equations:
        for expr, pos in equation.rhs.preorder_iter():
            if isinstance(expr, Times) and all(isinstance(operand, Symbol) or ((is_inverse(operand) or is_transpose(operand)) and isinstance(operand.operand, Symbol)) for operand in expr.operands):
                ops = set()
                for op in expr.operands:
                    if isinstance(op, Symbol):
                        ops.add(op)
                    else:
                        ops.add(op.operand)
                # We can't remove cases with one operand only because of (X^T X)^-1 X^T y, where QR does lead to a solution.
                if len(ops) > 1 and ops <= operands_to_factor:
                    yield set(op.name for op in ops)