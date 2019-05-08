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

        Returns the number of new nodes.
        """

        ###############
        # Apply kernels.

        new_nodes = []

        for node in self.active_nodes:

            transformed = self.TR_kernels(node.equations)

            for equations, edge_label, original_equations in transformed:
                equations = equations.remove_identities()

                new_nodes.extend(self.create_nodes(node, (equations, edge_label, original_equations)))

        return new_nodes

    def TR_kernels(self, equations):

        transformed_expressions = []

        for eqn_idx, equation in enumerate(equations):
            if not isinstance(equation.rhs, ae.Symbol):
                for eqns_variant in generate_variants(equations, eqn_idx):
                    pos, op_type = process_next_simple(eqns_variant[eqn_idx])

                    if op_type == OperationType.times and is_explicit_inversion(eqns_variant[eqn_idx][pos]):
                        transformed_expressions.extend(self.TR_matrix_chain(eqns_variant, eqn_idx, pos, True))
                    elif op_type == OperationType.unary:
                        transformed_expressions.extend(self.TR_unary_kernels(eqns_variant, eqn_idx, pos))
                    else:
                        te_reductions = self.TR_reductions(eqns_variant, eqn_idx, (1,))
                        transformed_expressions.extend(te_reductions)
                        if not te_reductions and not find_operands_to_factor(equations, eqn_idx):
                            # only use unary kernels if nothing else can be done
                            transformed_expressions.extend(self.TR_unary_kernels(eqns_variant, eqn_idx, (1,)))

                break

        return transformed_expressions

    def TR_reductions(self, equations, eqn_idx, initial_pos):

        transformed_expressions = []

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

                transformed_expressions.append((equations_copy, (matched_kernel,), equations))

        return transformed_expressions
