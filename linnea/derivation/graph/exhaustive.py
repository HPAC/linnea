from ...algebra.expression import Operator, Symbol
from ...algebra.transformations import simplify
from ...algebra.representations import to_SOP
from ...algebra.validity import check_validity
from ...algebra.properties import Property as properties

from . import base
from .utils import generate_variants, \
                   process_next, process_next_simple, \
                   OperationType, ExpressionType, \
                   is_explicit_inversion, \
                   DS_step, \
                   is_blocked

from .. import special_properties # TODO why is this necessary here?

import operator
import matchpy
import copy
import math

from ... import temporaries
from ... import config
from ..utils import select_optimal_match

collections_module = config.import_collections()

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

            new_nodes = self.DS_tricks()
            # TODO could this be done better with logging?
            self.print_DS_numbered("Nodes added (tricks):", new_nodes, self.level_counter)
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

        # if self.verbose
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
            if DS_step.kernels in node.applied_DS_steps:
                continue
            else:
                node.add_applied_step(DS_step.kernels)

            transformed = self.TR_kernels(node.equations)

            for equations, edge_label, original_equations in transformed:
                equations = equations.remove_identities()

                new_nodes.extend(self.create_nodes(node, (equations, edge_label, original_equations)))

            # Active node stops being active.
            # If no transformations were applied, it's a dead end.
            inactive_nodes.append(node)
  
        for node in inactive_nodes:
            self.active_nodes.remove(node)

        self.add_active_nodes(new_nodes)

        return len(new_nodes)

    def TR_kernels(self, equations):

        for eqn_idx, equation in enumerate(equations):
            if not isinstance(equation.rhs, Symbol):
                for eqns_variant in generate_variants(equations, eqn_idx):
                    pos, op_type = process_next_simple(eqns_variant[eqn_idx])

                    if op_type == OperationType.times and is_explicit_inversion(eqns_variant[eqn_idx][pos]):
                        for matrix_chain in self.TR_matrix_chain(eqns_variant, eqn_idx, pos, True):
                            yield matrix_chain
                    else:
                        reduction_yielded = False
                        for reduction in self.TR_reductions(eqns_variant, eqn_idx, (1,)):
                            reduction_yielded = True
                            yield reduction
                        if not reduction_yielded:
                            for unary_kernel in self.TR_unary_kernels(eqns_variant, eqn_idx, (1,)):
                                yield unary_kernel

                break

    def TR_reductions(self, equations, eqn_idx, initial_pos):

        initial_node = equations[eqn_idx][initial_pos]

        for node, _pos in initial_node.preorder_iter():
            pos = initial_pos + _pos

            for grouped_kernels in collections_module.reduction_MA.match(node).grouped():

                kernel, substitution = select_optimal_match(grouped_kernels)
                
                matched_kernel = kernel.set_match(substitution, True, CSE_rules=False)
                if is_blocked(matched_kernel.operation.rhs):
                    continue

                evaled_repl = matched_kernel.replacement

                new_equation = matchpy.replace(equations[eqn_idx], pos, evaled_repl)

                equations_copy = equations.set(eqn_idx, new_equation)
                equations_copy = equations_copy.to_normalform()

                temporaries.set_equivalent(equations[eqn_idx].rhs, equations_copy[eqn_idx].rhs)

                yield (equations_copy, (matched_kernel,), equations)
