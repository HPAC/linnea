from .base import derivation

from .utils import generate_variants

import math

class PropertyGraph(derivation.DerivationGraphBase):

    def derivation(self):
        number_of_nodes = 0
        while number_of_nodes != len(self.nodes):
            # print(number_of_nodes, len(self.nodes))
            number_of_nodes = len(self.nodes)

            self.DS_tricks()
            self.DS_collapse_nodes()
            self.DS_factorizations()
            self.DS_collapse_nodes()

            # self.to_dot_file()

    def DS_factorizations(self):
        """Applies factorizations.

        This function is specifically desinged for applying factorization when
        applying transformations to expression with special properties.

        Important:
        - The metric is set to (infinity, zero) to make sure that
          nothing is ever discarded
        - It is assumed that there is only one single equation.

        """
        new_nodes = []
        inactive_nodes = []

        for node in self.active_nodes:

            transformed = []
            for eqns_variant in generate_variants(node.equations, 0):
                transformed.extend(self.TR_factorizations(eqns_variant, 0, [1], [math.inf, 0]))

            new_nodes.extend(self.create_nodes(node, *transformed))

            # Active node stops being active.
            # If transformations were applied, it's actually not a leaf anymore,
            # i.e. it has successors.
            # If no transformations were applied, it's a dead end.
            inactive_nodes.append(node)
        
        for node in inactive_nodes:
            self.active_nodes.remove(node)

        self.add_active_nodes(new_nodes)
        # self.active_nodes.extend(new_nodes)

        return len(new_nodes)