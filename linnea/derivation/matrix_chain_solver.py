from ..algebra.expression import Times

from . import blocked_operations
from .. import config
collections_module = config.import_collections()

from .graph import matrix_chain

import math

import matchpy

from . import utils as dgu

class MatrixChainNotComputable(Exception):
    pass

class MatrixChainSolver(object):
    """docstring for MatrixChainSolver"""
    def __init__(self, expr):
        self.expr = expr
        self.n = len(self.expr.operands) # TODO remove?

        n = len(self.expr.operands)

        self.costs = [[math.inf for _ in range(n)] for _ in range(n)]
        self.sol = [[None for _ in range(n)] for _ in range(n)]
        self.tmps = [[None for _ in range(n)] for _ in range(n)]
        self.matched_kernels = [[None for _ in range(n)] for _ in range(n)]

        for i in range(n):
            self.costs[i][i] = 0
            self.tmps[i][i] = self.expr.operands[i]

        self._solve()
        self.matched_kernels = list(self._constuct_solution(0, n-1))
        self.tmp = self.tmps[0][n-1]

    def _solve(self):
        """Solves the matrix chain problem.

        This algorithm also takes indexed operands into account. Requires that the
        sizes of operands are known.

        For an explanation of the algorithm, have a look at "Introduction to
        Algorithms" by Cormen, Leiserson and Rivest.
        """

        n = len(self.expr.operands)

        for length in range(1, n): # length of the product
            for i in range(n - length):
                j = i + length

                solution_found = False
                solution_kernel = None
                solution_match = None

                for k in range(i, j):

                    if not self.tmps[i][k] or not self.tmps[k+1][j]:
                        # If one of those two temporaries does not exist, the
                        # corresponding sub-chain was not computable. Thus, we can
                        # not proceed here.
                        # This does not mean that the entire chain is not computable.
                        continue

                    product = Times(self.tmps[i][k], self.tmps[k+1][j])

                    if blocked_operations.is_blocked(product):
                        # print("blocked", product)
                        continue


                    kernel, match = dgu.select_optimal_match(collections_module.matrix_chain_DN.match(product))
                    if not kernel:
                        continue

                    # print(kernel.signature, kernel.pattern)

                    # Here, the indices of the temporary are computed without
                    # computing the temporary because
                    # - computing the temporary is fairly expensive
                    # - if cost > cost[i][j], the temporary is not used anyway
                    idx_range = product.index_range()

                    single_cost = kernel.cost(match)*idx_range
                    cost = self.costs[i][k] + self.costs[k+1][j] + single_cost
                    if (cost < self.costs[i][j]):
                        solution_kernel = kernel
                        solution_match = match
                        self.costs[i][j] = cost
                        self.sol[i][j] = k
                        solution_found = True

                if solution_found:
                    # create_tmp is called outside of the k loop. This way, the
                    # number of calls to create_tmp can be reduced from n^3 in the
                    # worst case to at most n^2. Notice that in some cases, this
                    # has absolutely no effect because n^3 is not reached anyway.

                    # Important: tmp.equiv_expr is not set in general and not used
                    # for deriving properties because it significanlty affects
                    # performance (there are calls to to_canonical and symplify).
                    # Instead, it is only set for the temporary that represents the
                    # entire expression.

                    # It looks a lot like this does not really affect performance anymore.

                    # If we heavily rely on equivalent expressions, they always
                    # have to be set here.
                    if i == 0 and j == n-1:
                        matched_kernel = solution_kernel.set_match(solution_match, False, False, False, True, self.expr)
                        # print(solution_match, self.expr)
                    else:
                        matched_kernel = solution_kernel.set_match(solution_match, False, False, False, False, None)
                    # matched_kernel = solution_kernel.set_match(solution_match, False)

                    self.tmps[i][j] = matched_kernel.replacement
                    self.matched_kernels[i][j] = matched_kernel

        if not self.tmps[0][n-1]:
            # If there is no temporary for the entire chain, then no solution was
            # found.
            # I don't think it's possible to detect ealier that no solution will be
            # found. Even if some parts can not be computed at all (i.e. some
            # tmps[i][j]) remain empty, it's still possible that the entire chain can
            # be solved.
            raise MatrixChainNotComputable("no solution found")

    def _constuct_solution(self, i, j):
        """ Yields the matched kernels necessary to compute this matrix chain.

        """
        if i != j:
            yield from self._constuct_solution(i, self.sol[i][j])
            yield from self._constuct_solution(self.sol[i][j]+1, j)
            yield self.matched_kernels[i][j]
