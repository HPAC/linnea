from ...algebra import expression as ae

from . import base
from .utils import generate_variants, \
                   process_next_simple, \
                   OperationType, \
                   is_explicit_inversion, \
                   find_operands_to_factor, \
                   is_blocked, \
                   DS_step


from .. import special_properties # TODO why is this necessary here?

import matchpy

from ... import temporaries
from ... import config
from ..utils import select_optimal_match

collections_module = config.import_collections()

class DerivationGraph(base.derivation.DerivationGraphBase):

    def DS_kernels(self):
        """applies all kernels to all active nodes and creates new nodes


        """

        new_nodes = []
        for node in self.active_nodes:
            for equations, edge_label, original_equations in self.TR_kernels(node.equations):
                new_nodes.extend(self.create_nodes(node, (equations, edge_label, original_equations)))

        return new_nodes

    def DFS_kernels_constructive(self, node, equations):
        for equations, edge_label, original_equations in self.TR_kernels_DFS_constructive(equations):
            yield self.create_node(node, equations, edge_label, original_equations, previous_DS_step=DS_step.kernels)

    def TR_kernels_DFS_constructive(self, equations):
        equation, eqn_idx = equations.process_next()
        if not equation:
            return
        pos, op_type = process_next_simple(equation)

        if op_type == OperationType.times:
            yield from self.TR_matrix_chain(equations, eqn_idx, pos, is_explicit_inversion(equation[pos]))
        elif op_type == OperationType.plus:
            yield from self.TR_addition(equations, eqn_idx, pos)
        # elif op_type == OperationType.none or (op_type == OperationType.none and not find_operands_to_factor(equations, eqn_idx)):
        #     # only use unary kernels if nothing else can be done
        #     yield from self.TR_unary_kernels(equations, eqn_idx, pos)

    def DFS_kernels(self, node, equations):
        for equations, edge_label, original_equations in self.TR_kernels_DFS(equations):
            yield self.create_node(node, equations, edge_label, original_equations, previous_DS_step=DS_step.kernels)


    def TR_kernels(self, equations):

        equation, eqn_idx = equations.process_next()
        if not equation:
            return
        for eqns_variant in generate_variants(equations, eqn_idx):
            pos, op_type = process_next_simple(eqns_variant[eqn_idx])

            if op_type == OperationType.times and is_explicit_inversion(eqns_variant[eqn_idx][pos]):
                yield from self.TR_matrix_chain(eqns_variant, eqn_idx, pos, True)
            else:
                # yield_from can't be used in this case, because we need to know if a reduction was yielded
                reduction_yielded = False
                for reduction in self.TR_reductions(eqns_variant, eqn_idx, (1,)):
                    reduction_yielded = True
                    yield reduction
                if not reduction_yielded and not find_operands_to_factor(equations, eqn_idx):
                    # only use unary kernels if nothing else can be done
                    yield from self.TR_unary_kernels(eqns_variant, eqn_idx, (1,))


    def TR_kernels_DFS(self, equations):

        equation, eqn_idx = equations.process_next()
        if not equation:
            return
        pos, op_type = process_next_simple(equation)

        if op_type == OperationType.times and is_explicit_inversion(equation[pos]):
            yield from self.TR_matrix_chain(equations, eqn_idx, pos, True)
        else:
            # yield_from can't be used in this case, because we need to know if a reduction was yielded
            reduction_yielded = False
            for reduction in self.TR_reductions(equations, eqn_idx, (1,)):
                reduction_yielded = True
                yield reduction
            if not reduction_yielded:
                # only use unary kernels if nothing else can be done
                yield from self.TR_unary_kernels(equations, eqn_idx, (1,))


    def TR_reductions(self, equations, eqn_idx, initial_pos):

        initial_node = equations[eqn_idx][initial_pos]

        for node, _pos in initial_node.preorder_iter():
            pos = initial_pos + _pos

            for grouped_kernels in collections_module.reduction_MA.match(node).grouped():

                kernel, substitution = select_optimal_match(grouped_kernels)
                
                matched_kernel = kernel.set_match(substitution, True)
                if is_blocked(matched_kernel.operation.rhs):
                    continue

                evaled_repl = matched_kernel.replacement

                new_equation = matchpy.replace(equations[eqn_idx], pos, evaled_repl)

                equations_copy = equations.set(eqn_idx, new_equation)
                equations_copy = equations_copy.to_normalform()

                temporaries.set_equivalent(equations[eqn_idx].rhs, equations_copy[eqn_idx].rhs)
                # remove_identities has to be called after set_equivalent because
                # after remove_identities, eqn_idx may not be correct anymore
                equations_copy = equations_copy.remove_identities()

                yield (equations_copy, (matched_kernel,), equations)
