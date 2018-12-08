from . import expression as ae

from .properties import Property as properties

from . import transformations as at

import itertools
import copy
import operator

def to_SOP(expr):
    """Brings an expression to the "sum of products" form.

    TODO think about implementing the following
    Note that this function does not consider scalar sums. Scalar sums are
    treated as single scalars. The reason is that there is probably no case
    where transforming (alpha + beta)*A into alpha*A + beta*A is useful.

    IMPORTANT: This function modifies expr. Since it is possible that the
    outermost operator changes, it is possible that the variable originally
    passed to this function points to a wrong node in the expression tree.
    To avoid any problems, use this fuction as follows:

    expr = to_SOP(expr)
    """

    # TODO should this be part of simplify()?

    # if isinstance(expr, Times) and any(isinstance(operand, Plus) and not operand.has_property(properties.SCALAR) for operand in expr.operands):
    if isinstance(expr, ae.Times) and any(isinstance(operand, ae.Plus) for operand in expr.operands):
        factors = []
        for operand in expr.operands:
            # REMARK: It's possible to do
            # operand = to_SOP(operand)
            # here instead of in the return line.
            # I would say that conceptually, it's more elegant, because this
            # version starts transforming to SOP from the inside out (i.e.
            # bottom to top in the expression tree), while the other version
            # works in the opposite direction. Intuitively, I would say that
            # the bottom-up approach should be faster because the expressions
            # that are copied are simpler, but it's not the case.

            # if isinstance(operand, Plus) and not operand.has_property(properties.SCALAR):
            if isinstance(operand, ae.Plus):
                factors.append(operand.operands)
            else:
                factors.append([operand])

        # Products have to be copied to ensure that the expression tree is
        # actually a tree. If common subexpression in the expression tree are
        # references to the same object, bad things happen when modifying one
        # occurence of that expression (all other occurences are modified as
        # well).
        return ae.Plus(*(to_SOP(ae.Times(*product)) for product in itertools.product(*factors)))
        # return Plus(factors)

    elif isinstance(expr, ae.Operator):
        return type(expr)(*map(to_SOP, expr.operands))
    else:
        return expr


def _sequence_intersection(sequences):
    """Computes all intersections of all sequences.

    sequences has to be a list of sorted lists. This functions computes the
    intersection of all those lists.

    The idea is to simultaneously traverse all lists, collecting elements that
    appear in all lists. This is done by maintaining a list of pointers,
    where each pointer points to one element of a different list (the current
    elements). Initially, all pointers point to the first element.

    1. The smallest element of the current elements is determined.
    2. If all elements are equal to the smallest element, this element is part
       of the intersection.
    3. The pointers pointing to those elements are moved forward by one.
    4. This process is repeated, starting at 1., until all lists are completely
       traversed.

    The intersection is returned as a list of lists. Each inner list is one set
    of indices of one element of the intersection.

    Example:
    Input: [["a", "b", "c", "d"],["b", "c", "d"],["a", "c", "d"]]
    Output: [[2, 1, 1], [3, 2, 2]]

    Clearly, the intersection consists of two elements len(output) == 2.
    The first element is at position 2 in the first sequence, 1 in the second
    sequence and 1 in the third sequence [2, 1, 1].
    """

    n = len(sequences)

    pointer = [0 for _ in range(n)]
    not_end = [True for _ in range(n)]

    intersections = []

    while True:

        # print("P:", pointer)

        # Finding the index k of the first sequences that has some unprocessed
        # elements left.
        k = not_end.index(True)
        min_val = sequences[k][pointer[k]]
        min_sequences = [k]

        # It is sufficient to start at k+1 because for all i < k, not_end[i] is
        # False.
        for i in range(k+1, n):
            if not_end[i]:
                val = sequences[i][pointer[i]]
                # print(val, min_val)
                if val < min_val:
                    # print("smaller")
                    min_sequences = [i]
                    min_val = val
                elif val == min_val:
                    # print("equal")
                    min_sequences.append(i)

        # print(min_sequences)

        if len(min_sequences) == n:
            intersections.append(tuple(copy.deepcopy(pointer)))

        # Advancing pointers pointing to minimal elements
        for i in min_sequences:
            if pointer[i] + 1 == len(sequences[i]):
                not_end[i] = False
            else:
                pointer[i] += 1

        # print(not_end)
        if not any(not_end):
            break
        
        # print(min_val)
        # print(min_sequences)

    return intersections

def to_POS(expr, side="l"):
    """Tries to transform expr into a product of sums (POS).

    For linear algebra expressions, the product of sums is not unique. This
    function tries to transforms the given expression into a product of sums
    as follows. For a sums of products, it factors out all expression on one
    side of those product until that is no longer possible. Then, it does the
    same on the with expression on the other side of the remaining products.
    Which side is used first is specified by the optional argument side, which
    has to be either "l" (left side is considered first, default) or "r" (right
    side first). Note that those to options may lead to different expressions
    and that the result may not be completely in POS form.

    For a sum of products, if a non-scalar is factored out, all scalars that
    show up in all those products are factored out as well.

    """

    if side == "l":
        other_side = "r"
    elif side == "r":
        other_side = "l"
    else:
        print("Error in to_POS.")

    expr = _to_POS(expr, side)
    expr = _to_POS(expr, other_side)
    return expr

def _keyfunc_left(term):
    """keyfunc for _to_POS

    Used for factoring out the left-most common factor of a multiple products.
    """
    return term[2][0]

def _keyfunc_right(term):
    """keyfunc for _to_POS

    Used for factoring out the right-most common factor of a multiple products.
    """
    return term[2][-1]

def _to_POS(expr, side):
    """Tries to transform expr into a product of sums (POS).

    This function does the actual factoring. The idea is to recursively traverse
    the expression tree and, when a sum is encountered, factor out one
    non-scalar expression (potentially plus some additional scalars). Then, the
    function is recursively called on the remaining sum.

    Side (either "l" or "r") determines whether a non-scalar on the left side or
    on the right side is factored out.

    Notice that a POS is for linear algebra expressions is not unique. To get
    all possible forms, it would be necessary to call _to_POS recursively both
    with "l" and "r", becaue each one could, in some cases, lead to a different
    expression.
    """
    # Should this function also be called on the non-scalar that is factored
    # out?

    # Only use this function on non-scalar expressions. On scalar expressions,
    # the keyfuncs don't work (because they try to access one element of the
    # list containing non-scalar expressions).
    
    # TODO does this make it faster?
    if not isinstance(expr, ae.Plus) and not isinstance(expr, ae.Symbol):
        return type(expr)(*[_to_POS(operand, side) for operand in expr.operands])
    elif expr.has_property(properties.SCALAR) or isinstance(expr, ae.Symbol):
        # This function doesn't do anything with scalar
        return expr

    if side == "l":
        idx = 0
        keyfunc = _keyfunc_left
    elif side == "r":
        idx = -1
        keyfunc = _keyfunc_right
    else:
        print("Error in _to_POS.")
    
    new_operands = []

    terms = []
    # Here, a term is a operand of a sum. terms contains tuples
    # of three elements. If the term is Times, the tuple contains:
    # - the expression itself (Times)
    # - a list of scalar operands
    # - a list of non-scalar operands
    # Otherwise, it contains:
    # - the expression itself
    # - an empty list
    # - a list containing the expression (thus, it is a list of non-scalar
    #   operands. This is one of the reasons why this function should only be
    #   called on non-scalar expressions.)
    for operand in expr.operands:
        if isinstance(operand, ae.Times):
            scalars, non_scalars = operand.split_operands()
            terms.append((operand, scalars, non_scalars))
        else:
            terms.append((operand, [], [operand]))

    # Terms are grouped either by the first or the last element in the list of 
    # non-scalar operands.
    terms.sort(key=keyfunc)
    for key, group in itertools.groupby(terms, key=keyfunc):
        group = list(group)
        # print(key, group)
        if len(group) > 1:
            # In this case, there are multiple terms that have a factor in common.
            common_factor = [key]
            scalar_sequences = []
            non_scalar_sequences = []
            for _, scalars, non_scalars in group:
                scalar_sequences.append(scalars)
                del non_scalars[idx]
                non_scalar_sequences.append(non_scalars)    
            # scalar_sequences = [scalars for _, scalars, _ in group]

            # Here, we check if those terms that have a non-scalar factor in common
            # also have scalars in common. If yes, they are factored out as well.
            if all(len(scalars) > 0 for scalars in scalar_sequences):
                map(operator.methodcaller("sort"), scalar_sequences)
                intersection = _sequence_intersection(scalar_sequences)
                
                if intersection:
                    scalar_common_factors = [scalar_sequences[0][positions[0]] for positions in intersection]

                    # Following from the way the intersections are constructed,
                    # positions with smaller indices always come first (i.e. 
                    # the positions are ordered).
                    # If we deleted in the orginial order, after the first
                    # deletion the following indices would be wrong.
                    for positions in reversed(intersection):
                        for n, scalars in enumerate(scalar_sequences):
                            del scalars[positions[n]]

                    scalar_common_factors.sort()
                    common_factor = scalar_common_factors + common_factor
            
            # Determining if the identity matrix has to be inserted. It has to
            # be inserted if there is a summand in the rest that is not a
            # scalar:
            #     alpha*A*B + beta*B = (alpha*A + beta*I)*B
            # It does not have to be inserted if the entire rest is scalar:
            #     alpha*B + beta*B = (alpha + beta)*B
            use_identity = any(non_scalar_sequences)

            if use_identity:
                for expr_, scalars, non_scalars in group:
                    if non_scalars:
                        identity_size = (non_scalars[0].rows, non_scalars[-1].columns)

            # Constructing the new terms with common factors removed.         
            rest = []
            for expr_, scalars, non_scalars in group:
                scalars.sort()
                if len(non_scalars) == 0:
                    if scalars:
                        if use_identity:
                            # New rest is: Times(alpha, beta, I)
                            operands = scalars + [ae.IdentityMatrix(*identity_size)]
                            rest.append(ae.Times(*operands))
                        else:
                            if len(scalars) == 1:
                                # New rest is: alpha
                                rest.append(scalars[0])
                            else:
                                # New rest is: Times(alpha, beta)
                                operands = scalars
                                rest.append(ae.Times(*operands))
                    else:
                        if use_identity:
                            # New rest is: I
                            rest.append(ae.IdentityMatrix(*identity_size))
                        else:
                            # New rest is: 1
                            rest.append(ae.ConstantScalar(1))
                else:
                    if scalars:
                        # New rest is: Times(alpha, beta, A, B)
                        operands = scalars + non_scalars
                        rest.append(ae.Times(*operands))
                    else:
                        if len(non_scalars) == 1:
                            # New rest is: A
                            rest.append(non_scalars[0])
                        else:
                            # New rest is: Times(A, B)
                            operands = non_scalars
                            rest.append(ae.Times(*operands))

            rest_expr = _to_POS(ae.Plus(*rest), side)
            if side == "l":
                new_operands.append(at.simplify(ae.Times(*(common_factor + [rest_expr]))))
            elif side == "r":
                new_operands.append(at.simplify(ae.Times(*([rest_expr] + common_factor))))
        else:
            # In this case, there is no common non-scalar factor, so the
            # original expression is still a operand of this plus.
            new_operands.append(group[0][0])
    if len(new_operands) == 1:
        return new_operands[0]
    else:
        return type(expr)(*new_operands)
