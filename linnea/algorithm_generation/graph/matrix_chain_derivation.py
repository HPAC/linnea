"""Matrix Chain Derivation Graph

This module was used for the CGO 2018 paper to do experitments only with the
matrix chain algorithm.
"""

from ...algebra import expression as ae

from ...import config

from .. import reductions

from .base import equations

import math


class DerivationGraph(equations.EquationsGraph):

    def derivation(self):

        old_verbosity = config.verbosity
        config.set_verbosity(0)

        self.best_first_search(time_limit=math.inf)

        config.set_verbosity(old_verbosity)


    def successor_generator(self, node):
        yield from self.DFS_matrix_chain(node, node.equations)


    def DFS_matrix_chain(self, node, equations):
        for new_equations, edge_label in self.TR_matrix_chain(equations):
            yield self.create_node(node, new_equations, edge_label, equations)


    def TR_matrix_chain(self, equations):
        if isinstance(equations[0].rhs, ae.Times):
            yield from reductions.apply_matrix_chain_algorithm(equations, 0, (1,))
        else:
            yield from reductions.apply_unary_kernels(equations, 0, (1,))
