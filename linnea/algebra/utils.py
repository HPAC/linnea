
import itertools

def _inverse_bandwidth(bandwidth, size):
    """Computes the bandwidth of the inverse of a matrix.

    """
    lb, ub = bandwidth
    n, m = size
    if n == m:
        if lb > 0:
            lb = n-1
        if ub > 0:
            ub = m-1
        return (lb, ub)
    else:
        if lb == 0 and ub == 0:
            # The pseudoinverse of a rectangular diagonal matrix is still
            # diagonal (the elements on the inverse are inverted)
            return bandwidth
        else:
            # Otherwise, nothing can be said about the result
            return (n-1, m-1)

# @profile
def scalar_subexpressions(expr):
    """Returns a set of all scalar subexpression in a product.

    This function returns a list of pairs of integers. Each pair represents one
    scalar subexpression, where the first integer is the position of the
    beginning of the scalar expression and the second integer the position of
    the end.

    Example: x^T A x A       -> {(0, 2)}
             x^T A x x^T x A -> {(0, 4)}
             x^T alpha x     -> {(0, 2), (1, 1)}

    """
    # We first find pairs of factors where the first
    # one has one row and the second one has one column (think of finding
    # matching parenthesis).
    #
    # Example: x^T A x A results in scalar_expressions = {(0, 2)}
    #          0 and 2 denote where a scalar expressios begins and ends

    scalar_expressions = set()
    stack = []
    factors = expr.operands # for legibility

    for n, factor in enumerate(factors):
        size = factor.size
        if size[0] == 1:
            stack.append(n)
        if size[1] == 1 and stack:
            scalar_expressions.add((stack.pop(), n))

    # To deal with consecutive scalar expressions, for example
    # x^T A x x^T x A, we have to fuse consecutive pairs.
    # Example: {(0, 2), (3, 4)} becomes {(0, 4)}
    # Remark: This is similar to (but not the same as) computing the
    # transitive closure.
    # TODO is it possible to do this more elegantly by sorting?
    # TODO Is it even necessary to do that?

    while True:
        new_pair = ()
        delete1 = ()
        delete2 = ()

        # it's not possible to fuse more than two pairs per iteration
        for pair1, pair2 in itertools.product(scalar_expressions, scalar_expressions):
            if pair1[1]+1 == pair2[0]:
                new_pair = (pair1[0], pair2[1])
                delete1 = pair1
                delete2 = pair2
                break

        if new_pair==():
            break

        scalar_expressions.add(new_pair)
        scalar_expressions.remove(delete1)
        scalar_expressions.remove(delete2)

    return scalar_expressions