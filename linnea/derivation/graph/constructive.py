from ...algebra import expression as ae
from ...algebra import transformations as at
from ...algebra.representations import to_SOP
from ...algebra.validity import check_validity
from ...algebra.properties import Property as properties

from ... import config

from . import base
from .utils import generate_variants, \
                   process_next, process_next_simple, \
                   OperationType, ExpressionType, \
                   is_explicit_inversion, \
                   DS_step, find_operands_to_factor

from .. import special_properties # TODO why is this necessary here?

import operator
import matchpy
import math

class DerivationGraph(base.derivation.DerivationGraphBase):

    def derivation(self, solution_nodes_limit=math.inf, iteration_limit=100, merging=True, dead_ends=True):

        check_validity(self.root.equations)
        self.root.equations = self.root.equations.to_normalform()
        self.root.equations.infer_lhs_properties()

        self.init_temporaries(self.root.equations)

        self.root.metric = self.root.equations.metric()

        # for testing and debugging
        trace = []

        terminal_nodes = []
        new_nodes_per_iteration = 0

        for i in range(iteration_limit):

            new_nodes_per_iteration = 0
            all_new_nodes = []


            new_nodes = self.DS_tricks()
            all_new_nodes.extend(new_nodes)
            # TODO could this be done better with logging?
            self.print_DS_numbered("Nodes added (tricks):", len(new_nodes), self.level_counter)
            trace.append(len(new_nodes))
            new_nodes_per_iteration += len(new_nodes)

            new_nodes = self.DS_factorizations()
            all_new_nodes.extend(new_nodes)
            self.print_DS_numbered("Nodes added (fact):", len(new_nodes), self.level_counter)
            trace.append(len(new_nodes))
            new_nodes_per_iteration += len(new_nodes)         

            new_nodes = self.DS_CSE_replacement()
            all_new_nodes.extend(new_nodes)
            self.print_DS_numbered("Nodes added (CSE):", len(new_nodes), self.level_counter)
            trace.append(len(new_nodes))
            new_nodes_per_iteration += len(new_nodes)         

            new_nodes = self.DS_kernels()
            all_new_nodes.extend(new_nodes)
            self.print_DS_numbered("Nodes added (kernels):", len(new_nodes), self.level_counter)
            trace.append(len(new_nodes))
            new_nodes_per_iteration += len(new_nodes)

            self.active_nodes = all_new_nodes

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

            terminal_nodes = self.terminal_nodes()

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

        self.print("{:-<34}".format(""))
        self.print_DS("Solution nodes:", len(terminal_nodes))
        self.print_DS("Number of nodes:", len(self.nodes))

        data = self.root.equations.get_data()
        self.print("Data: {}".format(data))
        if terminal_nodes:
            _, cost = self.shortest_path()
            self.print("Best solution: {:.3g}".format(cost))
            self.print("Intensity: {:.3g}".format(cost/data))

        # from ... import temporaries
        # print("\n".join(["{}: {}".format(k, v) for k, v in temporaries._equivalent_expressions.items()]))
        # print("######")
        # print("\n".join(["{}: {}".format(k, v) for k, v in temporaries._table_of_temporaries.items()]))
        
        #print(self.equivalence_rules)
        # produce output

        return trace

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

        transformed_expressions = []

        for eqn_idx, equation in enumerate(equations):
            if not isinstance(equation.rhs, ae.Symbol):
                for eqns_variant in generate_variants(equations, eqn_idx):
                    pos, op_type = process_next_simple(eqns_variant[eqn_idx])

                    if op_type == OperationType.times:
                        expr = eqns_variant[eqn_idx][pos]
                        transformed_expressions.extend(self.TR_matrix_chain(eqns_variant, eqn_idx, pos, is_explicit_inversion(expr)))
                    elif op_type == OperationType.plus:
                        transformed_expressions.extend(self.TR_addition(eqns_variant, eqn_idx, pos))
                    elif op_type == OperationType.unary or (op_type == OperationType.none and not find_operands_to_factor(equations, eqn_idx)):
                        # only use unary kernels if nothing else can be done
                        transformed_expressions.extend(self.TR_unary_kernels(eqns_variant, eqn_idx, pos))

                break

        return transformed_expressions