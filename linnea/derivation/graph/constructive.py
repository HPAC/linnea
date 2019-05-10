from ...algebra import expression as ae

from . import base
from .utils import generate_variants, \
                   process_next_simple, \
                   OperationType, \
                   is_explicit_inversion, \
                   find_operands_to_factor

from .. import special_properties # TODO why is this necessary here?

class DerivationGraph(base.derivation.DerivationGraphBase):

    def DS_kernels(self):
        """applies all kernels to all active nodes and creates new nodes

        Returns the number of new nodes.
        """

        new_nodes = []

        for node in self.active_nodes:

            transformed = self.TR_kernels(node.equations)

            for equations, edge_label, original_equations in transformed:
                equations = equations.remove_identities()

                new_nodes.extend(self.create_nodes(node, (equations, edge_label, original_equations)))

        return new_nodes


    def TR_kernels(self, equations):

        for eqn_idx, equation in enumerate(equations):
            if not isinstance(equation.rhs, ae.Symbol):
                for eqns_variant in generate_variants(equations, eqn_idx):
                    pos, op_type = process_next_simple(eqns_variant[eqn_idx])

                    if op_type == OperationType.times:
                        expr = eqns_variant[eqn_idx][pos]
                        yield from self.TR_matrix_chain(eqns_variant, eqn_idx, pos, is_explicit_inversion(expr))
                    elif op_type == OperationType.plus:
                        yield from self.TR_addition(eqns_variant, eqn_idx, pos)
                    elif op_type == OperationType.none or (op_type == OperationType.none and not find_operands_to_factor(equations, eqn_idx)):
                        # only use unary kernels if nothing else can be done
                        yield from self.TR_unary_kernels(eqns_variant, eqn_idx, pos)
                break
