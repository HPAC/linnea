from ..algebra.expression import Times

from . import blocked_operations
from .. import config
collections_module = config.import_collections()

from .graph import matrix_chain

import math
import itertools

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
        # TODO self.matched_kernels shouldn't be two entirely different things
        self.matched_kernels = list(itertools.chain.from_iterable(self._constuct_solution(0, n-1)))
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

                DN_solution_found = False
                DN_solution_kernel = None
                DN_solution_match = None

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

                    DN_solution = False
                    kernel, match = dgu.select_optimal_match(collections_module.matrix_chain_DN.match(product))
                    if kernel:
                        _cost = kernel.cost(match)
                        DN_solution = True
                    else:
                        # The matrix chain graph is only used if product can
                        # con be computed with a single kernel.
                        mc_graph = matrix_chain.MatrixChainGraph(product)
                        mc_graph.derivation()
                        matched_kernels, _cost, tmp = mc_graph.optimal_algorithm()
                        # mc_graph.to_dot_file(name="matrix_chain_graph.gv")
                        if not matched_kernels:
                            # mc_graph.to_dot_file(name="matrix_chain_graph.gv")
                            continue

                    cost = self.costs[i][k] + self.costs[k+1][j] + _cost
                    if (cost < self.costs[i][j]):
                        if DN_solution:
                            # DN_solution_kernel = kernel
                            # DN_solution_match = match
                            # DN_solution_found = True

                            matched_kernel = kernel.set_match(match, False)

                            self.tmps[i][j] = matched_kernel.replacement
                            self.matched_kernels[i][j] = [matched_kernel]
                        else:
                            self.matched_kernels[i][j] = matched_kernels
                            self.tmps[i][j] = tmp

                        self.costs[i][j] = cost
                        self.sol[i][j] = k

                # if DN_solution_found:
                #     # Calling set_match here is a performance optimization. It
                #     # is a fairly expensive function, and before this point,
                #     # only the cost is needed, which can be obtained separately,
                #     # and much cheaper.
                #     # This way, it is only called if the results are really
                #     # needed. The number of calls is reduced from n^3 in the
                #     # worst case to at most n^2. However, it is possible that
                #     # this has absolutely no effect because n^3 is not reached
                #     # anyway.
                #     # This optimization is not possible when using the matrix
                #     # chain graph because in that case, it's not possible to
                #     # obtain the cost without doing all of the other
                #     # computations.
                #     matched_kernel = DN_solution_kernel.set_match(DN_solution_match, False)

                #     self.tmps[i][j] = matched_kernel.replacement
                #     self.matched_kernels[i][j] = [matched_kernel]

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
