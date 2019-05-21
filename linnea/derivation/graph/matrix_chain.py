from ... import config

from ..utils import select_optimal_match

from .base import expression as egb

import matchpy
import operator

collections_module = config.import_collections()

class MatrixChainGraph(egb.ExpressionGraphBase):

    def derivation(self):
        nodes = [self.root]
        while nodes:
            new_nodes = []
            for node in nodes:
                new_nodes.extend(self.create_nodes(node, *self.TR_matrix_chain_kernels(node.expression)))
                new_nodes.extend(self.create_nodes(node, *self.TR_unary_kernels(node.expression)))
            nodes = new_nodes

            self.DS_merge_nodes()


    def DS_merge_nodes(self):
        """Merges redundant nodes in the derivation graph.

        Returns the number of removed nodes.

        """
        hashtable = dict()
        for node in self.nodes:
            hashtable.setdefault(node.get_payload(), []).append(node)

        remove = []
        for group in hashtable.values():
            remaining_node = group.pop()
            for node in group:
                remaining_node.merge(node)
            remove.extend(group)

        self.remove_nodes(remove)

        return len(remove)


    def TR_matrix_chain_kernels(self, expression):
        kernel, substitution = select_optimal_match(collections_module.matrix_chain_DN.match(expression))

        if kernel:
            matched_kernel = kernel.set_match(substitution, False)
            transformed_expression = matched_kernel.replacement

            return [(transformed_expression, (matched_kernel,))]
        return []


    def TR_unary_kernels(self, expression):
        transformed_expressions = []

        # iterate over all subexpressions
        for node, pos in expression.preorder_iter():
            kernel, substitution = select_optimal_match(collections_module.unary_kernel_DN.match(node))

            if kernel:
                matched_kernel = kernel.set_match(substitution, False)
                transformed_expression = matchpy.replace(expression, pos, matched_kernel.replacement)
                transformed_expressions.append((transformed_expression, (matched_kernel,)))

        return transformed_expressions
