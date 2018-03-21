from .base import expression as egb


from ... import config

collections_module = config.import_collections()

import matchpy
import operator

from ..utils import select_optimal_match

class MatrixChainGraph(egb.ExpressionGraphBase):

    def derivation(self):
        _break = False
        while self.active_nodes:
            # self.DS_matrix_chain_kernels()
            # self.DS_unary_kernels()
            # # stop if there is a solution node? Does it ever pay off?

            new_nodes = []
            inactive_nodes = []

            for node in self.active_nodes:
                transformed = self.TR_matrix_chain_kernels(node.expression)

                new_nodes.extend(self.create_nodes(node, *transformed))
                # if any(map(operator.methodcaller("is_terminal"), new_nodes)):
                #     _break = True
                #     break

                transformed = self.TR_unary_kernels(node.expression)

                new_nodes.extend(self.create_nodes(node, *transformed))

                inactive_nodes.append(node)

            if _break:
                _break = False
                break

            for node in inactive_nodes:
                self.active_nodes.remove(node)

            self.add_active_nodes(new_nodes)

            self.DS_merge_nodes()

    # TODO remove?
    def DS_matrix_chain(self):
        new_nodes = []
        inactive_nodes = []

        for node in self.active_nodes:

            transformed = self.TR_matrix_chain_kernels(node.expression)
            transformed.extend(self.TR_unary_kernels(node.expression))

            new_nodes.extend(self.create_nodes(node, *transformed))

            # Active node stops being active.
            # If transformations were applied, it's actually not a leaf anymore,
            # i.e. it has successors.
            # If no transformations were applied, it's a dead end.
            inactive_nodes.append(node)
        
        for node in inactive_nodes:
            self.active_nodes.remove(node)

        self.add_active_nodes(new_nodes)

        return len(new_nodes)


    def TR_matrix_chain_kernels(self, expression):
        kernel, substitution = select_optimal_match(collections_module.matrix_chain_DN.match(expression))

        if kernel:
            
            # print("match", kernel.pattern, substitution)
                
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
                # print([(kernel.pattern, substitution) for kernel, substitution in matches])
                
                matched_kernel = kernel.set_match(substitution, False)

                evaled_repl = matched_kernel.replacement
                transformed_expression = matchpy.replace(expression, pos, evaled_repl)

                transformed_expressions.append((transformed_expression, (matched_kernel,)))

        return transformed_expressions