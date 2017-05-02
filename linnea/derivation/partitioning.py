from ..algebra.expression import Plus, Times, Symbol, Matrix, \
                                      BlockedExpression, Equal, \
                                      Inverse, Conjugate, Transpose, \
                                      ConjugateTranspose, \
                                      InverseTranspose, InverseConjugate, \
                                      InverseConjugateTranspose

from ..algebra.properties import Property as properties

from ..utils import clamp

import itertools
import functools
import copy

TOP = dict() # table of partitionings

# Template for the name of partitioned operands.
block_name = "{0}[{1}{2}]"

def clear():
    TOP.clear()

def window(seq, n=2):
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    result = tuple(itertools.islice(it, n))
    if len(result) == n:
        yield result    
    for elem in it:
        result = result[1:] + (elem,)
        yield result

def apply_partitioning(expr):
    # This has to be checked here because apply_partitioning can and should not
    # be called on the parts.
    if isinstance(expr, BlockedExpression):
        return expr
    if not isinstance(expr, Symbol):
        expr.operands = [apply_partitioning(operand) for operand in expr.operands]
    if isinstance(expr, Symbol):
        if expr.partitioning == (set(), set()):
            return expr

        rows, cols = expr.partitioning
        nrows = len(rows) + 1
        ncols = len(cols) + 1

        row_IDs = []
        if nrows == 1:
            row_IDs = [""]
        elif nrows == 2:
            row_IDs = ["T", "B"]
        elif nrows == 3:
            row_IDs = ["T", "M", "B"]
        else:
            row_IDs = [str(i) for i in range(nrows)]

        col_IDs = []
        if ncols == 1:
            col_IDs = [""]
        elif ncols == 2:
            col_IDs = ["L", "R"]
        elif ncols == 3:
            col_IDs = ["L", "M", "R"]
        else:
            col_IDs = [str(i) for i in range(ncols)]

        symbol_rows, symbol_cols = expr.size

        if nrows == 1:
            row_sizes = [symbol_rows]
        else:
            rows = sorted(rows)
            row_sizes = [rows[0]]
            for n1, n2 in window(rows):
                row_sizes.append(n2-n1)
            row_sizes.append(symbol_rows - rows[-1])

        if ncols == 1:
            col_sizes = [symbol_cols]
        else:
            cols = sorted(cols)
            col_sizes = [cols[0]]
            for n1, n2 in window(cols):
                col_sizes.append(n2-n1)
            col_sizes.append(symbol_cols - cols[-1])

        lb, ub = expr.bandwidth
        # print(row_IDs)
        # print(col_IDs)
        # print(row_sizes)
        # print(col_sizes)
        row_partitioning = list(itertools.accumulate(row_sizes))
        row_partitioning.insert(0, 0)
        col_partitioning = list(itertools.accumulate(col_sizes))
        col_partitioning.insert(0, 0)
        # print(row_partitioning)
        # print(col_partitioning)

        operands = []
        for i in range(nrows):
            row = []
            for j in range(ncols):
                # print("block", i, j)
                _lb = lb + col_partitioning[j] - row_partitioning[i]
                _ub = ub + row_partitioning[i] - col_partitioning[j]
                # print(_lb, -col_sizes[j], row_sizes[i]-1)
                _lb = clamp(_lb, -col_sizes[j], row_sizes[i]-1)
                _ub = clamp(_ub, -row_sizes[i], col_sizes[j]-1)
                size = (row_sizes[i], col_sizes[j])
                # If the number of diagonals is zero, generate zero matrix
                if _lb + _ub + 1 <= 0:
                    block = Zero(size)
                else:
                    if expr.has_property(properties.IDENTITY):
                        block = Identity(size)
                    else:
                        name = block_name.format(expr.name, row_IDs[i], col_IDs[j])
                        block = Matrix(name, size, expr.indices)
                        block.bandwidth = (_lb, _ub)
                # print(_lb, _ub)
                row.append(block)
            operands.append(row)
        return BlockedExpression(operands)
    elif isinstance(expr, Times):
        if any(isinstance(operand, BlockedExpression) for operand in expr.operands):

            scalars, non_scalars = expr.split_operands()

            # Multiply
            factors = []
            for operand in non_scalars:
                if isinstance(operand, BlockedExpression):
                    factors.append(operand.operands)
                else:
                    factors.append([[operand]])

            product_operands = functools.reduce(multiply_blocked_matrices, factors)

            if scalars:
                product_operands = [[Times(block, *scalars) for block in row] for row in product_operands]

            # Test if result is just a single block. If yes, just return that
            # single block, not BlockedExpression.
            if len(product_operands)==1 and len(product_operands[0])==1:
                return product_operands[0][0]
            else:
                return BlockedExpression(product_operands)
        else:
            return expr

    elif isinstance(expr, Plus):
        if isinstance(expr.operands[0], BlockedExpression):
            operands = expr.operands
            # print(operands)
            # TODO it seriously bothers me that this works. operands is a list
            # of BlockedExpressions, they are iterable. However, for Transpose,
            # map_thread doesn't work, even though that's how it's implemented
            # in BlockeExpression_utils.py (I suspect some of that stuff does
            # not work.)
            # print(operands)
            return BlockedExpression( map_thread( Plus, operands, 2 ))
        else:
            return expr
    elif isinstance(expr, Equal):
        if isinstance(expr.operands[0], BlockedExpression):
            operands = expr.operands
            # TODO same here
            return BlockedExpression( map_thread( Equal, operands, 2 ))
        else:
            return expr
    elif isinstance(expr, Transpose):
        if isinstance(expr.operand, BlockedExpression):
            blocked_expr = expr.operand
            # list(zip(*[...])) transposes the list
            # since zip returns tuples, not lists, map(list, ...) is added
            blocked_expr.operands = list(map(list, zip(*[[Transpose(elem) for elem in row] for row in blocked_expr.operands])))
            return blocked_expr
        else:
            return expr
    else:
        # print("apply_partitioning can not be applied to ", expr)
        return expr
    # TODO Further unary operators missing

def multiply_blocked_matrices( m1, m2 ):
    m2_transposed = list(zip(*m2))
    product = []
    if len(m2) == 1:
        # Number of rows of matrix 2 is one (that is the "inner" dimension of
        # the product). In this case, there is no Plus. (For example, this could
        # be an outer product. This also covers the case where m1 and/or m2 is not
        # blocked, that is, m1 and/or m2 is a single matrix.)
        for row1 in m1:
            row = []
            for col2 in m2_transposed:
                row.append(Times(copy.deepcopy(row1[0]), copy.deepcopy(col2[0])))
            product.append(row)
    elif len(m1) == 1 and len(m2_transposed) == 1:
        # m1 has one row and m2 has one column, so this is an inner product. In
        # this case, nothing has to be copied since each element is just used
        # once.
        product = [[Plus(*(Times(elem1, elem2) for elem1, elem2 in zip(m1[0], m2_transposed[0])))]]
    else:
        for row1 in m1:
            row = []
            for col2 in m2_transposed:
                row.append(Plus(*(Times(copy.deepcopy(elem1), copy.deepcopy(elem2)) for elem1, elem2 in zip(row1, col2))))
            product.append(row)
    return product

def propagate_partitioning(expr):
    """Propagates all partitionings of operands in expr.

    """
    while _propagate_partitioning(expr):
        pass

def _propagate_partitioning(expr):
    """One iteration of propagating partitionings in expr.

    Returns True if something changed, False if not.
    """
    ret = False
    if isinstance(expr, Plus):
        scalars, non_scalars = expr.split_operands()
        ret = ret or combine_sets(expr.partitioning[0], *[operand.partitioning[0] for operand in non_scalars])
        ret = ret or combine_sets(expr.partitioning[1], *[operand.partitioning[1] for operand in non_scalars])
    elif isinstance(expr, Times):
        scalars, non_scalars = expr.split_operands()
        if non_scalars:
            ret = ret or combine_sets(expr.partitioning[0], non_scalars[0].partitioning[0])
            ret = ret or combine_sets(expr.partitioning[1], non_scalars[-1].partitioning[1])
            for c1, c2 in window(non_scalars):
                ret = ret or combine_sets(c1.partitioning[1], c2.partitioning[0])
    elif isinstance(expr, Transpose) or isinstance(expr, ConjugateTranspose):
        ret = ret or combine_sets(expr.partitioning[0], expr.operand.partitioning[1])
        ret = ret or combine_sets(expr.partitioning[1], expr.operand.partitioning[0])
    elif isinstance(expr, Conjugate):
        ret = ret or combine_sets(expr.partitioning[0], expr.operand.partitioning[0])
        ret = ret or combine_sets(expr.partitioning[1], expr.operand.partitioning[1])
    # TODO missing: Inverse, InverseConjugate, InverseTranspose, InverseConjugateTranspose
    # inverse with anything but (set(), set()) should probably cause an error
    elif isinstance(expr, Symbol):
        if any(expr.has_property(prop) for prop in [properties.SYMMETRIC, properties.TRIANGULAR, properties.DIAGONAL]): # what if diagonal (...) and not square?
            ret = ret or combine_sets(*expr.partitioning)
        if not expr.name in TOP:
            TOP[expr.name] = (set(), set())
        partitioning = TOP[expr.name]
        ret = ret or combine_sets(partitioning[0], expr.partitioning[0])
        ret = ret or combine_sets(partitioning[1], expr.partitioning[1])
    elif isinstance(expr, BlockedExpression):
        if expr.partitioning == (set(), set()):
            row_partitioning = set(itertools.accumulate([row[0].size[0] for row in expr.operands][:-1]))
            col_partitioning = set(itertools.accumulate([block.size[1] for block in expr.operands[0]][:-1]))
            expr.partitioning = (row_partitioning, col_partitioning)
            ret = True

    elif isinstance(expr, Equal):
        ret = ret or combine_sets(expr.lhs.partitioning[0], expr.rhs.partitioning[0])
        ret = ret or combine_sets(expr.lhs.partitioning[1], expr.rhs.partitioning[1])
    else:
        # print("missing operator in _propagate_partitioning:", expr)
        return False

    if not isinstance(expr, BlockedExpression) and not isinstance(expr, Symbol):
        # Using any() is not possible here because it short-circuits, but
        # _propagate_partitioning has to be called on all operands.
        for operand in expr.operands:
            ret = ret or _propagate_partitioning(operand)
    return ret

def combine_sets(*sets):
    """Updates all sets with the union of all sets.

    This function takes an arbitray number of sets. It updates all sets such
    that each set is the union of all original sets. Sets are updated
    in place, that is, the sets passed to this function are modified.

    The function returns True if any set is modified and False is nothing
    changed.
    """
    u = set.union(*sets)
    # An alternative way to test if something changes would be to test if all
    # sets are equal. If this is the case, nothing changes. However, all sets
    # have to be updated anyway, so they have to be touched anyway. len is
    # probably really cheap, so it probably doesn't matter.
    # Remember that the sets are usually tiny (0-2 elements).
    change = False
    for s in sets:
        l1 = len(s)
        s.update(u)
        l2 = len(s)
        if l1 != l2:
            change = True
    return change
