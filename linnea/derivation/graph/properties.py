from ... import config

from .base import derivation
from .utils import generate_variants

import math
import itertools

class PropertyGraph(derivation.DerivationGraphBase):

    def derivation(self):

        old_verbosity = config.verbosity
        config.set_verbosity(0)

        self.best_first_search(time_limit=math.inf)

        config.set_verbosity(old_verbosity)
        
    def successor_generator(self, node):

        funs = [self.DFS_factorizations, self.DFS_tricks]
        generators = [fun(node, eqns) for eqns, fun in itertools.product(generate_variants(node.equations), funs)]
        yield from itertools.chain.from_iterable(generators)