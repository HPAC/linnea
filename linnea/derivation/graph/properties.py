from ... import config

from .base import derivation
from .utils import generate_variants

import itertools

class PropertyGraph(derivation.DerivationGraphBase):

    def derivation(self):

        old_verbosity = config.verbosity
        config.set_verbosity(0)

        """
        In most cases, this search finishes in much less than a second. However,
        there might be (rare) cases where it takes much longer. To avoid
        unexpected behavior of the derivation, we limit the search to 1 second.
        This may prevent special_properties from working properly, i.e. when
        factorizations or tricks are applied,
        - merging may work less often, and
        - some properties are not preserved.
        """
        self.best_first_search(time_limit=1)

        config.set_verbosity(old_verbosity)
        
    def successor_generator(self, node):

        funs = [self.DFS_factorizations, self.DFS_tricks]
        generators = [fun(node, eqns) for eqns, fun in itertools.product(generate_variants(node.equations), funs)]
        yield from itertools.chain.from_iterable(generators)