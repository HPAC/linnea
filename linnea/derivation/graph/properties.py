from .base import derivation

from .utils import generate_variants

import math

class PropertyGraph(derivation.DerivationGraphBase):

    def derivation(self):
        number_of_nodes = 0
        while number_of_nodes != len(self.nodes):
            number_of_nodes = len(self.nodes)

            self.DS_tricks()
            self.DS_merge_nodes()
            self.DS_factorizations()
            self.DS_merge_nodes()