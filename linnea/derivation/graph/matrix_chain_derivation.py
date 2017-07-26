from ...algebra import expression as ae
from ...algebra import transformations as at

# from ...algebra.expression import Operator, Symbol
# from ...algebra.transformations import simplify
from ...algebra.representations import to_SOP
from ...algebra.consistency import check_consistency
from ...algebra.validity import check_validity
from ...algebra.properties import Property as properties

from . import base
from .utils import generate_variants, \
                   process_next, process_next_simple, \
                   OperationType, ExpressionType

from .. import special_properties

import operator
import matchpy

class DerivationGraph(base.derivation.DerivationGraphBase):

    def derivation(self, max_algos=1, max_iterations=5, verbose=True):
        # TODO default values?

        self.verbose = verbose

        # init_partitiong and delete_partitioning is only done for performance
        # reasons. The additional attribute "partitioning" makes copying a lot
        # more expensive.
        # TODO this should by unnecessary by now

        # for equation in self.root.equations:
        #     equation.init_partitioning()

        # self.root.equations.apply_partitioning()

        # for equation in self.root.equations:
        #     equation.delete_partitioning()

        # change order with partitioning?
        # self.root.equations.resolve_dependencies()
        self.root.equations.replace_auxiliaries()

        check_validity(self.root.equations)
        for eqn_idx, equation in enumerate(self.root.equations):
            self.root.equations[eqn_idx] = to_SOP(at.simplify(equation))
            for node, pos in equation.rhs.preorder_iter():
                # I think it's not necessary to store properties both for
                # expr and Inverse(expr) if expr is SPD. Think about that. 
                if not (isinstance(node, ae.Operator) and node.arity is matchpy.Arity.unary):
                    if node.has_property(properties.SPD):
                        special_properties.add_expression(node, [properties.SPD])

            # TODO check_consistency for Equations?
            check_consistency(equation)

        self.root.metric = self.root.equations.metric()

        suc_deriv = 0

        # for testing and debugging
        trace = []


        new_nodes = self.DS_kernels()
        self.print_DS_numbered("Nodes added (kernels):", new_nodes, self.level_counter)
        trace.append(new_nodes)


        # TODO order of merge and prune?
        merged_nodes = self.DS_collapse_nodes()
        self.print_DS("Nodes merged:", merged_nodes)
        trace.append(merged_nodes)


        terminal_nodes = list(filter(operator.methodcaller("is_terminal"), self.nodes))

        if len(terminal_nodes) >= 0:
            self.print("Specified number of algorithms found.")
            return []
        elif not self.active_nodes:
            self.print("No further derivations possible.")
            return []

            # self.to_dot_file("counter")
            # print("Leaves", [node.id for node in self.active_nodes])
            # print("Nodes", [node.id for node in self.nodes])


        self.print("{:-<30}".format(""))
        self.print_DS("Solution nodes:", len(terminal_nodes))
        self.print_DS("Number of nodes:", len(self.nodes))

        if terminal_nodes:
            self.print("Best solution: {:.3g}".format(min(map(operator.attrgetter("accumulated_cost"), terminal_nodes))))

        # if self.verbose:
        #     from ... import temporaries
        #     print("\n".join(["{}: {}".format(k, v) for k, v in temporaries._equivalent_expressions.items()]))
        #     print("\n".join(["{}: {}".format(k, v) for k, v in temporaries._table_of_temporaries.items()]))
        
        #print(self.equivalence_rules)
        # produce output

        return trace

    def DS_kernels(self):
        """applies all kernels to all active nodes and creates new nodes

        Returns the number of new nodes.
        """

        ###############
        # Apply kernels.

        new_nodes = []
        inactive_nodes = []

        for node in self.active_nodes:

            transformed = self.TR_kernels(node.equations, node.metric)

            for equations, metric, edge_label in transformed:

                new_nodes.extend(self.create_nodes(node, (equations, metric, edge_label)))

            # Active node stops being active.
            # If no transformations were applied, it's a dead end.
            inactive_nodes.append(node)
  
        for node in inactive_nodes:
            self.active_nodes.remove(node)

        self.add_active_nodes(new_nodes)

        return len(new_nodes)


    def TR_kernels(self, equations, metric):

        transformed_expressions = []

        transformed_expressions.extend(self.TR_matrix_chain(equations, 0, [1], metric))

        return transformed_expressions