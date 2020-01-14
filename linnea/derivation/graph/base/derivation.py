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
        
        algorithms = list(self.k_best_algorithms(algorithms_limit, algorithm_name, no_duplicates, k_best, pruning_factor))

        self.print_result("Number of algorithms:", len(algorithms))

        if code or derivation or experiment_code:
            output_path = os.path.join(config.output_code_path, output_name)
            code_path = os.path.join(output_path, config.language.name, subdir_name)
            experiment_path = os.path.join(output_path, config.language.name, subdir_name_experiments)
            if not os.path.exists(output_path):
                os.makedirs(output_path)
            else:
                # Removing existing algorithm files
                if code:
                    cgu.remove_files(code_path)
                    cgu.remove_files(os.path.join(code_path, "derivation"))
                if experiment_code:
                    cgu.remove_files(experiment_path)
                    cgu.remove_files(os.path.join(experiment_path, "derivation"))

        if code or derivation:
            for n, algorithm in enumerate(algorithms):

                algorithm_file_name = "".join([algorithm.name, config.filename_extension])
                derivation_file_name = os.path.join("derivation", "".join([algorithm.name, ".txt"]))

                if code:
                    file_name = os.path.join(code_path, algorithm_file_name)
                    cgu.to_file(file_name, algorithm.code_as_function())
                    if derivation:
                        file_name = os.path.join(code_path, derivation_file_name)
                        cgu.to_file(file_name, algorithm.derivation())

                if experiment_code:
                    file_name = os.path.join(experiment_path, algorithm_file_name)
                    cgu.to_file(file_name, algorithm.code_as_function(experiment=True))
                    if derivation:
                        file_name = os.path.join(experiment_path, derivation_file_name)
                        cgu.to_file(file_name, algorithm.derivation())

        if experiment_code:
            generate_experiment_code(output_name, self.input, algorithm_name, [1, 24], k_best, len(algorithms))

        return algorithms


    def k_best_algorithms(self, k, algorithm_name="algorithm{}", no_duplicates=True, k_best=False, pruning_factor=1.0):
        """Generates the k best algorithms.

        Algorithms are generated in the order of increasing cost. Less than k
        algortihms are generated if fewer exist.

        Args:
            k (int): The number of algorithms to generated.
            algorithm_name (string, optional): Format string used as template
                for the algortihm name. Needs to contain exactly one occurrence
                of "{}", which is replaced with a number. Defaults to
                "algorithm{}".
            no_duplicates (bool, optional): If True, duplicate algorithms are
                removed. This may affect the performance if there is a large
                number of duplicate algorithms. Defaults to True.
            k_best (bool, optional): Flag for the "k best" experiment. If True,
                only algorithms with a cost smaller than "pruning_factor*cost of
                best algorithm" are generated (at most k).
            pruning_factor (float, optional): If k_best is True, only algorithms
                with a cost smaller than "pruning_factor*cost of best algorithm"
                are generated (at most k). Is ignored if k_best is False.
                Defaults to 1.0.

        Yields:
            Algorithm: The k best algorithms.
        """
        algorithms = []
        known_algorithms = set()
        number_of_algorithms = 0
        min_cost = self.shortest_path()[1]
        for path, cost in self.k_shortest_paths(math.inf):
            if k_best and cost > pruning_factor*min_cost:
                break

            algorithm = self.path_to_algorithm(path, cost, algorithm_name.format(number_of_algorithms))
            if no_duplicates:
                if not algorithm in known_algorithms:
                    known_algorithms.add(algorithm)
                    yield algorithm
                    number_of_algorithms += 1
            else:
                yield algorithm
                number_of_algorithms += 1

            if number_of_algorithms == k:
                break


    def optimal_algorithm_to_str(self):
        path, cost = self.shortest_path()
        if path:
            algorithm = self.path_to_algorithm(path, cost, "algorithm")
            return algorithm.code_as_function()
        else:
            return None


    def path_to_algorithm(self, path, cost, algorithm_name):
        kernels_and_equations = []
        current_node = self.root
        for idx in path:
            edge_label = current_node.edge_labels[idx]
            kernels_and_equations.append(current_node.original_equations[idx])
            kernels_and_equations.extend(edge_label.matched_kernels)
            current_node = current_node.successors[idx]

        return cgu.Algorithm(algorithm_name, self.input, current_node.equations, kernels_and_equations, cost)

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
