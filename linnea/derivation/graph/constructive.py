from ...algebra import expression as ae
from ...algebra import transformations as at

# from ...algebra.expression import Operator, Symbol
# from ...algebra.transformations import simplify
from ...algebra.representations import to_SOP
from ...algebra.consistency import check_consistency
from ...algebra.validity import check_validity
from ...algebra.properties import Property as properties

from ... import config

from . import base
from .utils import generate_variants, \
                   process_next, process_next_simple, \
                   OperationType, ExpressionType, \
                   is_explicit_inversion, \
                   DS_step

from .. import special_properties

import operator
import matchpy
import math

class DerivationGraph(base.derivation.DerivationGraphBase):

    # def prune_functionXX(self, nodes):
    #     return nodes

    def prune_function00(self, nodes):
        return nodes

    def prune_function01(self, nodes):
        nodes.sort(key=operator.attrgetter("accumulated_cost"))
        width = math.ceil(len(nodes)*0.5)
        return nodes[:width]

    def prune_function02(self, nodes):

        def _f(level):
            if level < 4:
                return 1
            else:
                return 0.3

        nodes.sort(key=operator.attrgetter("accumulated_cost"))
        width = math.ceil(len(nodes)*_f(self.level_counter))
        return nodes[:width]

    def derivation(self, solution_nodes_limit=math.inf, iteration_limit=100, merging=True, dead_ends=True):
        # TODO default values?

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
        # self.root.equations.replace_auxiliaries()

        check_validity(self.root.equations)
        self.root.equations = self.root.equations.to_normalform()
        for eqn_idx, equation in enumerate(self.root.equations):
            for node, pos in equation.rhs.preorder_iter():
                # I think it's not necessary to store properties both for
                # expr and Inverse(expr) if expr is SPD. Think about that. 
                if not (isinstance(node, ae.Operator) and node.arity is matchpy.Arity.unary):
                    if node.has_property(properties.SPD):
                        special_properties.add_expression(node, [properties.SPD])

            # TODO check_consistency for Equations?
            check_consistency(equation)

        self.root.metric = self.root.equations.metric()

        # for testing and debugging
        trace = []

        new_nodes_per_iteration = 0

        for i in range(iteration_limit):

            new_nodes_per_iteration = 0

            new_nodes = self.DS_tricks()
            # TODO could this be done better with logging?
            self.print_DS_numbered("Nodes added (tricks):", new_nodes, self.level_counter)
            trace.append(new_nodes)
            new_nodes_per_iteration += new_nodes

            if merging and new_nodes:
                merged_nodes = self.DS_merge_nodes()
                self.print_DS("Nodes merged:", merged_nodes)
                trace.append(merged_nodes)

            new_nodes = self.DS_CSE_replacement()
            self.print_DS_numbered("Nodes added (CSE):", new_nodes, self.level_counter)
            trace.append(new_nodes)
            new_nodes_per_iteration += new_nodes

            if merging and new_nodes:
                merged_nodes = self.DS_merge_nodes()
                self.print_DS("Nodes merged:", merged_nodes)
                trace.append(merged_nodes)

            new_nodes = self.DS_factorizations()
            self.print_DS_numbered("Nodes added (fact):", new_nodes, self.level_counter)
            trace.append(new_nodes)
            new_nodes_per_iteration += new_nodes

            if new_nodes:
                if dead_ends:
                    dead_nodes = self.DS_dead_ends()
                    self.print_DS("Dead nodes:", dead_nodes)
                    trace.append(dead_nodes)

                if merging:
                    merged_nodes = self.DS_merge_nodes()
                    self.print_DS("Nodes merged:", merged_nodes)
                    trace.append(merged_nodes)

            new_nodes = self.DS_CSE_replacement()
            self.print_DS_numbered("Nodes added (CSE):", new_nodes, self.level_counter)
            trace.append(new_nodes)
            new_nodes_per_iteration += new_nodes
            
            if merging and new_nodes:
                merged_nodes = self.DS_merge_nodes()
                self.print_DS("Nodes merged:", merged_nodes)
                trace.append(merged_nodes)

            new_nodes = self.DS_kernels()
            self.print_DS_numbered("Nodes added (kernels):", new_nodes, self.level_counter)
            trace.append(new_nodes)
            new_nodes_per_iteration += new_nodes

            # DEBUG TODO remove
            # for node in self.active_nodes:
            #     for equation in node.equations:
            #         check_consistency(equation)

            # for node in self.active_nodes:
            #     rhs = node.expression.rhs
            #     for node, pos in rhs.iterate_preorder_with_positions():
            #         print("###")
            #         print(node)
            #         print(node.properties)
            #         print(node.false_properties)
                # for prop in properties.all:
                #     has_property automatically stores properties on the node
                #    node.has_property(prop)

            if new_nodes:
                if dead_ends:
                    dead_nodes = self.DS_dead_ends()
                    self.print_DS("Dead nodes:", dead_nodes)
                    trace.append(dead_nodes)

                if merging:
                    merged_nodes = self.DS_merge_nodes()
                    self.print_DS("Nodes merged:", merged_nodes)
                    trace.append(merged_nodes)

            mins = self.metric_mins()
            #print(mins)
            pruned_nodes = self.DS_prune(mins)
            self.print_DS("Nodes pruned:", pruned_nodes)
            trace.append(pruned_nodes)

            terminal_nodes = list(filter(operator.methodcaller("is_terminal"), self.nodes))

            if len(terminal_nodes) >= solution_nodes_limit:
                self.print("Specified number of algorithms found.")
                break
            elif not self.active_nodes or new_nodes_per_iteration == 0:
                self.print("No further derivations possible.")
                break

            # self.to_dot_file("counter")
            # print("Leaves", [node.id for node in self.active_nodes])
            # print("Nodes", [node.id for node in self.nodes])
        else:
            self.print("Iteration limit reached.")

        self.print("{:-<30}".format(""))
        self.print_DS("Solution nodes:", len(terminal_nodes))
        self.print_DS("Number of nodes:", len(self.nodes))

        if terminal_nodes:
            self.print("Best solution: {:.3g}".format(min(map(operator.attrgetter("accumulated_cost"), terminal_nodes))))

        # from ... import temporaries
        # print("\n".join(["{}: {}".format(k, v) for k, v in temporaries._equivalent_expressions.items()]))
        # print("\n".join(["{}: {}".format(k, v) for k, v in temporaries._table_of_temporaries.items()]))
        
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
            if DS_step.kernels in node.applied_DS_steps:
                """Why exactly is this necessary? How is a state in the graph
                reached where this is necessary? This may have to do with
                leaving nodes active: As a results, it's possible to have
                a node and it's successor be active at the same time.
                """
                inactive_nodes.append(node)
                continue
            else:
                node.add_applied_step(DS_step.kernels)

            transformed = self.TR_kernels(node.equations, node.metric)

            for equations, metric, edge_label in transformed:
                equations = equations.remove_identities()
                # TODO move that into TR_kernels (in an elegant way)?
                # TODO in the best case, identity equations do not contribute to the metric
                metric = equations.metric()

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

        for eqn_idx, equation in enumerate(equations):
            if not isinstance(equation.rhs, ae.Symbol):
                for eqns_variant in generate_variants(equations, eqn_idx):
                    pos, op_type = process_next_simple(eqns_variant[eqn_idx])

                    if op_type == OperationType.times:
                        expr = eqns_variant[eqn_idx][pos]
                        transformed_expressions.extend(self.TR_matrix_chain(eqns_variant, eqn_idx, pos, metric, is_explicit_inversion(expr)))
                    elif op_type == OperationType.plus:
                        transformed_expressions.extend(self.TR_addition(eqns_variant, eqn_idx, pos, metric))
                    elif op_type == OperationType.none:
                        transformed_expressions.extend(self.TR_unary_kernels(eqns_variant, eqn_idx, pos, metric))

                break

        return transformed_expressions