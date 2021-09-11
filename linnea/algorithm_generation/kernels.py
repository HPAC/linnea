from ..algebra.expression import Equal, Plus, Times, Index, Symbol
from ..utils import number_of_entries
from ..code_generation.utils import MatchedKernel

from .. import temporaries
from .. import config

from .utils import select_optimal_match

from .graph.utils import process_next, OperationType, \
                         is_explicit_inversion, is_blocked

from . import matrix_chain_solver

import copy
import itertools
import math
import matchpy

collections_module = config.import_collections()

class MatrixSumNotComputable(Exception):
    pass


def apply_kernels_constructive(equations):
    equation, eqn_idx = equations.process_next()
    if not equation:
        return
    pos, op_type = process_next(equation)

    if op_type == OperationType.times:
        yield from apply_matrix_chain_algorithm(equations, eqn_idx, pos, is_explicit_inversion(equation[pos]))
    elif op_type == OperationType.plus:
        yield from apply_sum_algorithm(equations, eqn_idx, pos)
    elif op_type == OperationType.unary:
        yield from apply_unary_kernels(equations, eqn_idx, pos)


def apply_kernels_exhaustive(equations):
    equation, eqn_idx = equations.process_next()
    if not equation:
        return

    pos, op_type = process_next(equation)

    if op_type == OperationType.unary:
        # There is no need to do anything in this case here because
        # apply_kernels_constructive does it already.
        return
    else:
        # yield_from can't be used in this case, because we need to know if a kernel was yielded
        kernel_yielded = False
        for kernel in apply_kernels(equations, eqn_idx, (1,)):
            kernel_yielded = True
            yield kernel
        if not kernel_yielded:
            # only use unary kernels if nothing else can be done
            yield from apply_unary_kernels(equations, eqn_idx, (1,))


def apply_kernels(equations, eqn_idx, initial_pos):

    initial_node = equations[eqn_idx][initial_pos]

    for node, _pos in initial_node.preorder_iter():
        pos = initial_pos + _pos

        for grouped_kernels in collections_module.kernel_MA.match(node).grouped():

            kernel, substitution = select_optimal_match(grouped_kernels)
            
            matched_kernel = kernel.set_match(substitution, True)
            if is_blocked(matched_kernel.operation.rhs):
                continue

            evaled_repl = matched_kernel.replacement

            new_equation = matchpy.replace(equations[eqn_idx], pos, evaled_repl)

            equations_copy = equations.set(eqn_idx, new_equation)
            equations_copy = equations_copy.remove_explicit_transposition(eqn_idx)
            equations_copy = equations_copy.to_normalform()
            equations_copy = equations_copy.remove_identities()

            yield (equations_copy, (matched_kernel,))


def apply_matrix_chain_algorithm(equations, eqn_idx, initial_pos, explicit_inversion=False):

    try:
        msc = matrix_chain_solver.MatrixChainSolver(equations[eqn_idx][initial_pos], explicit_inversion)
    except matrix_chain_solver.MatrixChainNotComputable:
        return

    replacement = msc.tmp
    matched_kernels = msc.matched_kernels

    new_equation = matchpy.replace(equations[eqn_idx], initial_pos, replacement)
    equations_copy = equations.set(eqn_idx, new_equation)
    equations_copy = equations_copy.remove_explicit_transposition(eqn_idx)
    equations_copy = equations_copy.to_normalform()
    equations_copy = equations_copy.remove_identities()        

    yield (equations_copy, matched_kernels)


def apply_sum_algorithm(equations, eqn_idx, initial_pos):

    # Note: For addition, we decided to only use a binary kernel, no
    #       variadic addition.

    new_expr, matched_kernels = decompose_sum(equations[eqn_idx][initial_pos])

    new_equation = matchpy.replace(equations[eqn_idx], initial_pos, new_expr)
    equations_copy = equations.set(eqn_idx, new_equation)
    equations_copy = equations_copy.to_normalform().remove_identities()

    yield (equations_copy, matched_kernels)


def apply_unary_kernels(equations, eqn_idx, initial_pos):

    initial_node = equations[eqn_idx][initial_pos]

    # iterate over all subexpressions
    for node, _pos in initial_node.preorder_iter():
        pos = initial_pos + _pos

        kernel, substitution = select_optimal_match(collections_module.unary_kernel_DN.match(node))

        if kernel:
            matched_kernel = kernel.set_match(substitution, False)
            if is_blocked(matched_kernel.operation.rhs):
                continue
            evaled_repl = matched_kernel.replacement

            equations_copy = equations.set(eqn_idx, matchpy.replace(equations[eqn_idx], pos, evaled_repl))
            equations_copy = equations_copy.to_normalform().remove_identities()

            yield (equations_copy, (matched_kernel,))


def decompose_sum(expr):
    matched_kernels = []

    # Note: For addition, we decided to only use a binary kernel, no
    #       variadic addition.

    while isinstance(expr, Plus):
        expr, matched_kernel = apply_kernel_with_context(expr, collections_module.addition_kernel_MA)
        if matched_kernel:
            matched_kernels.append(matched_kernel)
            continue

        expr, matched_kernel = apply_kernel_anywhere(expr, collections_module.transposition_kernel_DN)
        if matched_kernel:
            matched_kernels.append(matched_kernel)
            continue

        expr, matched_kernel = apply_kernel_anywhere(expr, collections_module.scaling_kernel_DN)
        if matched_kernel:
            matched_kernels.append(matched_kernel)
            continue

        """
        When
        - expr has more than one operand and
        - the lists of kernels are exhausted without finding a match
        then expr can not be computed.
        """
        raise MatrixSumNotComputable("{} is not computable.".format(expr))

    return expr, matched_kernels


def apply_kernel_with_context(expr, many_to_one_matcher):
    """Applies the optimal kernel to an expressions.

    This function applies the optimal matching kernel of the many_to_one_matcher
    to expr and returns the replacement, as well as the matched kernel.

    The functions only searches for matches at the root of expr, and it expects
    many_to_one_matcher to use the patterns with context variables
    (Kernel.pattern_with_context).

    Args:
        expr (Expression): The expression where the kernel should be applied.
        many_to_one_matcher (ManyToOneMatcher): Matcher for kernel patterns with
            context. The label has to be the kernel.

    Returns:
        A tuple containing
        - the replacement (Expression)
        - the matched kernel (MatchedKernel)
        If no match was found, the replacement is the input expression and
        instead of the matched kernel, None is returned.
    """

    optimal_match = (None, None)
    min_cost = math.inf

    for kernel, substitution in many_to_one_matcher.match(expr):
        cost = kernel.cost(substitution)
        if cost < min_cost:
            optimal_match = (kernel, substitution)

    kernel, substitution = optimal_match
    if kernel:
        matched_kernel = kernel.set_match(substitution, True)
        expr = matched_kernel.replacement
        return expr, matched_kernel
    return (expr, None)


def apply_kernel_anywhere(expr, discrimination_net):
    """Applies one kernel to an expressions.

    This function applies the first matching kernel of the discrimination_net to
    expr and returns the replacement, as well as the matched kernel.

    The functions searches for matches anywhere in expr, and it expects
    discrimination_net to use the pattern without context variables
    (Kernel.pattern).

    Args:
        expr (Expression): The expression where the kernel should be applied.
        discrimination_net (DiscriminationNet): Discrimination net for kernel
        patterns without context. The final_label has to be the kernel.

    Returns:
        A tuple containing
        - the replacement (Expression)
        - the matched kernel (MatchedKernel)
        If no match was found, the replacement is the input expression and
        instead of the matched kernel, None is returned.
    """
    for node, pos in expr.preorder_iter():
        for kernel, substitution in discrimination_net.match(node):
            matched_kernel = kernel.set_match(substitution, False)
            expr = matchpy.replace(expr, pos, matched_kernel.replacement)
            return expr, matched_kernel
    return (expr, None)


def decompose_sum_old(expr):
    """Decomposes a sum with indexed operands to reduce FLOP count.

    This function returns a tuple of three elements:
    - The first element is a temporary representing the result of this sum (to
      be used as a replacement for the original sum.)
    - The second element is a list of the operations needed to compute this sum
    - The third argument is a list of the cost of those operations.
    """

    # Index sets are ordered by size (including the union of all sets).
    # Example: A + B{i} + C{j} + D{ik}
    # results in
    # [
    #  [{}],
    #  [{i}, {j}],
    #  [{i, k}],
    #  [{i, j, k}],
    # ]

    all_indices = set([frozenset(operand.indices) for operand in expr.operands])
    union_of_all_sets = frozenset(expr.indices)
    all_indices.add(union_of_all_sets)
    all_indices = sorted(all_indices, key=lambda x: len(x))
    all_indices = [list(g) for k, g in itertools.groupby(all_indices, key=lambda x: len(x))]

    # This list is used to store all operands that we still have to add.
    operands = copy.copy(expr.operands)

    # output_operations = []
    # output_costs = []

    operations = []

    # We go from small index sets to the union of all sets
    for indices_list in all_indices:

        min_idx_range = math.inf 
        min_indices = None
        # if True:
        # For index list with more than one element (i.e. [{i}, {j}]), we search
        # for the one with the smallest index range (if sizes are known).
        if len(indices_list[0]) != 0:
            # Testing if current index set is not empty
            # TODO do I need to test that?
            for indices in indices_list:
                idx_range = Index.range(indices)
                if idx_range < min_idx_range:
                    min_idx_range = idx_range
                    min_indices = indices
        # elif len(indices_list[0]) != 0:
        #     # If sizes are not known, we use an arbitrary one.
        #     min_indices = indices_list[0]

        for indices in indices_list:
            if indices == min_indices or min_indices == None:
                # If the current indices are the one with the smallest index
                # range or there is just one index, we add all operands that have
                # subsets of indices.
                ops = [operand for operand in filter(lambda x: x.indices.issubset(indices), operands)]
            else:
                # Otherwise, we just add operands with the same indices.
                ops = [operand for operand in filter(lambda x: x.indices == indices, operands)]
            if len(ops) > 1:
                # Only if there is more than one operand a temporary is generated.
                plus_expr = Plus(*ops)
                if indices == union_of_all_sets:
                    # For the final sum, equiv_expr is set.
                    tmp = temporaries.create_tmp(expr, True)
                else:
                    tmp = temporaries.create_tmp(plus_expr, False)
                # The cost is computed as
                # number of nonzero entries of result * number of operands * index range
                # Operands that come with an addition scalar multiplication are 
                # taken into account by adding the number of those operands to
                # the number of actual operands.
                scalar_multiplications = len([op for op in ops if isinstance(op, Times)])
                cost = number_of_entries(tmp) * (len(ops) - 1 + scalar_multiplications) * tmp.index_range()
                operations.append(Equal(tmp, plus_expr))
                # pseudo_code_fragment.add_line(Equal(tmp, plus_expr), "Matrix sum", cost)
                # output_operations.append(Equal(tmp, plus_expr))
                # output_costs.append(cost)
                for operand in ops:
                    operands.remove(operand)
                operands.append(tmp)
    
    # print(clak.kernels.collections.matrix_sum.pattern_with_context)
    # print(expr)
    # matches = list(patternmatcher.functional.match(expr, clak.kernels.collections.matrix_sum.pattern_with_context))
    return (operations[-1].lhs, [MatchedKernel()])



