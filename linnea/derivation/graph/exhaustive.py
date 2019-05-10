from ...algebra import expression as ae

from . import base
from .utils import generate_variants, \
                   process_next_simple, \
                   OperationType, \
                   is_explicit_inversion, \
                   find_operands_to_factor, \
                   is_blocked


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

    def TR_kernels(self, equations):

        for eqn_idx, equation in enumerate(equations):
            if not isinstance(equation.rhs, ae.Symbol):
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

                break

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
                equations_copy = equations_copy.to_normalform().remove_identities()

                temporaries.set_equivalent(equations[eqn_idx].rhs, equations_copy[eqn_idx].rhs)

                yield (equations_copy, (matched_kernel,), equations)
