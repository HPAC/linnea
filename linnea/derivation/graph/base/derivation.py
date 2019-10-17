from ....algebra.expression import Symbol
from ....code_generation import utils as cgu
from ....code_generation.experiments.utils import generate_experiment_code

from .... import config

from ... import tricks
from ... import CSEs
from ... import factorizations
from ... import reductions

from ..utils import generate_variants, find_operands_to_factor, \
                    DS_step, is_dead_end, PriorityStack, process_next_simple, \
                    OperationType, is_explicit_inversion

from . import base

import itertools
import os.path
import math
import time

class DerivationGraphBase(base.GraphBase):
    """

    Graph for equations.
    """

    def __init__(self, input):
        super().__init__()
        self.input = input
        self.root = DerivationGraphNode(input)
        self.nodes = [self.root]


    def best_first_search(self, time_limit=60, merging=True, dead_ends=True, pruning_factor=1.0):

        self.root.generator = self.successor_generator(self.root)

        p_stack = PriorityStack()
        p_stack.put(0, self.root)

        existing_nodes = dict()
        t_start = time.perf_counter()
        best_solution = math.inf
        trace_data = []
        terminal_nodes = []

        print_interval = 2
        next_print = print_interval
        pruned_nodes = dict()

        while not p_stack.empty():
            prio, node = p_stack.get()

            if node.accumulated_cost > best_solution * pruning_factor:
                node.labels.append("pruned")
                pruned_nodes[node.id] = (prio, node)
                continue

            if dead_ends and all(is_dead_end(equations, node.factored_operands) for equations in generate_variants(node.equations)):
                node.labels.append("dead")
                continue

            nodes_merged = False
            new_terminal_node = False

            try:
                new_node = next(node.generator)
            except StopIteration:
                pass
            else:
                # TODO is there a good way to avoid code duplication here?
                if merging:
                    try:
                        existing_prio, existing_node = existing_nodes[new_node.equations]
                    except KeyError:
                        existing_nodes[new_node.equations] = (0, new_node)
                        if new_node.is_terminal():
                            terminal_nodes.append(new_node)
                            new_terminal_node = True
                        else:
                            new_node.generator = self.successor_generator(new_node)
                            p_stack.put(0, new_node)
                    else:
                        existing_node.merge(new_node)
                        self.remove_node(new_node)
                        nodes_merged = True
                else:
                    if new_node.is_terminal():
                        terminal_nodes.append(new_node)
                        new_terminal_node = True
                    else:
                        new_node.generator = self.successor_generator(new_node)
                        p_stack.put(0, new_node)

                p_stack.put(prio+1, node)
            
            if new_terminal_node or nodes_merged:
                for terminal_node in terminal_nodes:
                    if terminal_node.accumulated_cost < best_solution:
                        best_solution = terminal_node.accumulated_cost
                        t_elapsed = time.perf_counter() - t_start
                        trace_data.append((t_elapsed, terminal_node.accumulated_cost))
                        self.print_result_sep("New solution:", "{:.3g}".format(best_solution))

            # TODO how to improve this?
            # put this into a function
            if nodes_merged:
                remove = []
                for p_prio, p_node in pruned_nodes.values():
                    if p_node.accumulated_cost <= best_solution * pruning_factor:
                        p_node.labels.remove("pruned") 
                        p_stack.put(p_prio, p_node)
                        remove.append(p_node.id)
                for id in remove:
                    del pruned_nodes[id]
           
            t_elapsed = time.perf_counter() - t_start
            if t_elapsed > next_print:
                self.print_result_sep("Nodes:", len(self.nodes))
                next_print += print_interval

            if t_elapsed > time_limit:
                self.print("Time limit reached.")
                break

        else:
            self.print("No further derivations possible.")

        return trace_data, terminal_nodes


    def DFS_kernels_constructive(self, node, equations):
        for new_equations, edge_label in self.TR_kernels_constructive(equations):
            yield self.create_node(node, new_equations, edge_label, equations, previous_DS_step=DS_step.kernels)


    def DFS_kernels(self, node, equations):
        for new_equations, edge_label in self.TR_kernels(equations):
            yield self.create_node(node, new_equations, edge_label, equations, previous_DS_step=DS_step.kernels)


    def DFS_tricks(self, node, equations):
        for new_equations, edge_label in tricks.apply_tricks(equations):
            yield self.create_node(node, new_equations, edge_label, equations, previous_DS_step=DS_step.tricks)


    def DFS_CSE_replacement(self, node, equations):
        for new_equations in CSEs.find_CSEs(equations):
            yield self.create_node(node, new_equations, (), equations, previous_DS_step=DS_step.CSE)


    def DFS_factorizations(self, node, equations):
        """Factorization derivation step.

        This derivation step applies factorizations to the active nodes. For a
        description of how exactly factorizations are applied, see the docstring
        of TR_factorizations().
        """

        # find all matrices that need to factored
        operands_to_factor = find_operands_to_factor(equations)

        operands_to_factor -= node.factored_operands
        if operands_to_factor:
            # construct dict {operand name: list of matched kernels for all valid factorizations}
            # It is constructed here because it can be reused. 
            factorization_dict = factorizations.construct_factorization_dict(operands_to_factor)

            for new_equations, edge_label in factorizations.apply_factorizations(equations, operands_to_factor, factorization_dict):
                yield self.create_node(node, new_equations, edge_label, equations, factored_operands=operands_to_factor, previous_DS_step=DS_step.factorizations)


    def TR_kernels_constructive(self, equations):
        equation, eqn_idx = equations.process_next()
        if not equation:
            return
        pos, op_type = process_next_simple(equation)

        if op_type == OperationType.times:
            yield from reductions.apply_matrix_chain_algorithm(equations, eqn_idx, pos, is_explicit_inversion(equation[pos]))
        elif op_type == OperationType.plus:
            yield from reductions.apply_sum_algorithm(equations, eqn_idx, pos)


    def TR_kernels(self, equations):
        equation, eqn_idx = equations.process_next()
        if not equation:
            return
        pos, op_type = process_next_simple(equation)

        if op_type == OperationType.times and is_explicit_inversion(equation[pos]):
            yield from reductions.apply_matrix_chain_algorithm(equations, eqn_idx, pos, True)
        else:
            # yield_from can't be used in this case, because we need to know if a reduction was yielded
            reduction_yielded = False
            for reduction in reductions.apply_reductions(equations, eqn_idx, (1,)):
                reduction_yielded = True
                yield reduction
            if not reduction_yielded:
                # only use unary kernels if nothing else can be done
                yield from reductions.apply_unary_kernels(equations, eqn_idx, (1,))


    def create_node(self, predecessor, equations, matched_kernels, original_equations, factored_operands=None, previous_DS_step=None):
        new_node = DerivationGraphNode(equations, factored_operands, previous_DS_step)
        new_node.level = self.step_counter
        predecessor.set_labeled_edge(new_node, base.EdgeLabel(*matched_kernels), original_equations)
        self.nodes.append(new_node)
        return new_node


    def create_nodes(self, predecessor, *description, factored_operands=None, previous_DS_step=None):
        new_nodes = []
        # print(description)
        # if description:

        _factored_operands = predecessor.factored_operands
        if factored_operands:
            _factored_operands = _factored_operands.union(factored_operands)

        for equations, matched_kernels, original_equations in description:
            new_node = self.create_node(predecessor, equations, matched_kernels, original_equations, _factored_operands.copy(), previous_DS_step)
            new_nodes.append(new_node)
        return new_nodes


    def write_output(
            self,
            code=True,
            derivation=False,
            output_name="tmp",
            experiment_code=False,
            k_best=False,
            algorithms_limit=1,
            pruning_factor=1.0,
            graph=False,
            graph_style=config.GraphStyle.full,
            subdir_name="generated",
            subdir_name_experiments="experiments",
            algorithm_name="algorithm{}",
            no_duplicates=False
        ):

        if not config.output_code_path:
            raise config.OutputPathNotSet("Unable to write output: output_code_path not set.")

        if graph:
            self.write_graph(output_name, graph_style)
        
        algorithms = []
        known_algorithms = set()
        number_of_algorithms = 0
        min_cost = self.shortest_path()[1]
        for path, cost in self.k_shortest_paths(math.inf):
            if k_best and cost > pruning_factor*min_cost:
                break

            algorithm = self.path_to_algorithm(path, cost)
            if no_duplicates:
                if not algorithm in known_algorithms:
                    known_algorithms.add(algorithm)
                    algorithms.append(algorithm)
                    number_of_algorithms += 1
            else:
                algorithms.append(algorithm)
                number_of_algorithms += 1

            if number_of_algorithms == algorithms_limit:
                break

        self.print_result("Number of algorithms:", number_of_algorithms)

        if code or derivation or experiment_code:
            directory_name = os.path.join(config.output_code_path, output_name)
            if not os.path.exists(directory_name):
                os.makedirs(directory_name)
            else:
                # Removing existing algorithm files
                if code:
                    algorithms_dir_name = os.path.join(directory_name, config.language.name, subdir_name)
                    cgu.remove_files(algorithms_dir_name)
                    derivation_dir_name = os.path.join(directory_name, config.language.name, subdir_name, "derivation")
                    cgu.remove_files(derivation_dir_name)
                if experiment_code:
                    algorithms_dir_name = os.path.join(directory_name, config.language.name, subdir_name_experiments)
                    cgu.remove_files(algorithms_dir_name)
                    algorithms_dir_name = os.path.join(directory_name, config.language.name, subdir_name_experiments, "derivation")
                    cgu.remove_files(algorithms_dir_name)

        if code or derivation:
            for n, algorithm in enumerate(algorithms):

                if code:
                    cgu.algorithm_to_file(output_name, subdir_name, algorithm_name.format(n), algorithm.code(), algorithm.experiment_input, algorithm.experiment_output)
                    if derivation:
                        cgu.derivation_to_file(output_name, subdir_name, algorithm_name.format(n), algorithm.derivation())

                if experiment_code:
                    cgu.algorithm_to_file(output_name, subdir_name_experiments, algorithm_name.format(n), algorithm.code(), algorithm.experiment_input, algorithm.experiment_output, experiment=True)
                    if derivation:
                        cgu.derivation_to_file(output_name, subdir_name_experiments, algorithm_name.format(n), algorithm.derivation())

        if experiment_code:
            generate_experiment_code(output_name, self.input, algorithm_name, [1, 24], k_best, number_of_algorithms)


        return algorithms

    def optimal_algorithm_to_str(self):
        matched_kernels, cost, final_equations = self.optimal_algorithm()
        if matched_kernels:
            algorithm = cgu.Algorithm(self.input, final_equations, matched_kernels, cost)
            code = algorithm.code()
            # code_lines = code.splitlines()
            # del code_lines[1:3] # remove "using Base.LinAlg..."
            # code = "\n".join(code_lines)
            function_name = "linnea_function{:X}".format(hash(self.input))
            return function_name, cgu.algorithm_to_str(function_name, code, algorithm.experiment_input, algorithm.experiment_output)
        else:
            return None, None


    def path_to_algorithm(self, path, cost):
        kernels_and_equations = []
        current_node = self.root
        for idx in path:
            edge_label = current_node.edge_labels[idx]
            kernels_and_equations.append(current_node.original_equations[idx])
            kernels_and_equations.extend(edge_label.matched_kernels)
            current_node = current_node.successors[idx]

        return cgu.Algorithm(self.input, current_node.equations, kernels_and_equations, cost)

class DerivationGraphNode(base.GraphNodeBase):

    _counter = 0

    def __init__(self, equations=None, factored_operands=None, previous_DS_step=None):
        super().__init__(factored_operands, previous_DS_step)

        # IDs for dot output
        self.id = DerivationGraphNode._counter
        self.name = "".join(["node", str(self.id)])
        DerivationGraphNode._counter +=1
        self.equations = equations
        self.generator = None
        
    def get_payload(self):
        return self.equations

    def is_terminal(self):
        """Returns true if the expression of self is fully computed."""
        return all([isinstance(equation.rhs, Symbol) for equation in self.equations])
