from . import base

from ..utils import generate_variants, find_operands_to_factor, \
                    find_occurrences, InverseType, group_occurrences, \
                    DS_step, find_explicit_symbol_inverse, is_inverse, \
                    find_blocking_products

from ....algebra.expression import Symbol, Times

from ....algebra.equations import Equations
from ....algebra.consistency import check_consistency
from ....algebra.validity import check_validity
from ....utils import powerset

import matchpy
import itertools
import os.path
import math

from ... import special_properties
from ... import tricks
from ... import CSEs
from ... import matrix_chain_solver as mcs
from ...matrix_sum import decompose_sum

from ...utils import select_optimal_match
from ....code_generation import utils as cgu
from ....code_generation.experiments import operand_generation, runner, reference_code
from .... import config
from .... import temporaries


collections_module = config.import_collections()

class DerivationGraphBase(base.GraphBase):
    """

    Graph for equations.
    """

    def __init__(self, input):
        super().__init__()
        self.input = input
        self.root = DerivationGraphNode(input)
        self.active_nodes = [self.root]
        self.nodes = [self.root]


    def derivation(self, solution_nodes_limit=math.inf, iteration_limit=100, merging=True, dead_ends=True):

        check_validity(self.root.equations)
        self.root.equations = self.root.equations.to_normalform()
        self.root.equations.infer_lhs_properties()

        self.init_temporaries(self.root.equations)

        self.root.metric = self.root.equations.metric()

        # for testing and debugging
        trace = []

        terminal_nodes = []
        self.step_counter = 1

        for i in range(iteration_limit):

            new_nodes_per_iteration = 0
            all_new_nodes = []

            new_nodes = self.DS_tricks()
            all_new_nodes.extend(new_nodes)
            # TODO could this be done better with logging?
            self.print_DS_numbered("Nodes added (tricks):", len(new_nodes), self.step_counter)
            trace.append(len(new_nodes))
            self.step_counter += 1

            new_nodes = self.DS_factorizations()
            all_new_nodes.extend(new_nodes)
            self.print_DS_numbered("Nodes added (fact):", len(new_nodes), self.step_counter)
            trace.append(len(new_nodes))
            self.step_counter += 1         

            new_nodes = self.DS_CSE_replacement()
            all_new_nodes.extend(new_nodes)
            self.print_DS_numbered("Nodes added (CSE):", len(new_nodes), self.step_counter)
            trace.append(len(new_nodes))
            self.step_counter += 1         

            new_nodes = self.DS_kernels()
            all_new_nodes.extend(new_nodes)
            self.print_DS_numbered("Nodes added (kernels):", len(new_nodes), self.step_counter)
            trace.append(len(new_nodes))
            self.step_counter += 1

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
            elif not self.active_nodes or not all_new_nodes:
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


    def create_nodes(self, predecessor, *description, factored_operands=None, previous_DS_step=None):
        new_nodes = []
        # print(description)
        # if description:

        _factored_operands = predecessor.factored_operands
        if factored_operands:
            _factored_operands = _factored_operands.union(factored_operands)

        for equations, matched_kernels, original_equations in description:
            new_node = DerivationGraphNode(equations, predecessor, _factored_operands.copy(), previous_DS_step)
            new_node.level = self.step_counter
            self.nodes.append(new_node)
            new_nodes.append(new_node)
            predecessor.set_labeled_edge(new_node, base.EdgeLabel(*matched_kernels), original_equations)
        return new_nodes


    def write_output(self, code=True, derivation=False, output_name="tmp", experiment_code=False, algorithms_limit=1, graph=False, graph_style=config.GraphStyle.full, subdir_name="generated", algorithm_name="algorithm{}"):

        if not config.output_code_path:
            raise config.OutputPathNotSet("Unable to write output: output_code_path not set.")

        if graph:
            self.write_graph(output_name, graph_style)
        
        paths = []
        min_cost = self.shortest_path()[1]
        for path, cost in self.k_shortest_paths(algorithms_limit):
            if cost <= 1.5*min_cost:
                paths.append((path, cost))
            else:
                break

        number_of_algorithms = len(paths)
        if config.verbosity >= 1:
            self.print("Number of algorithms: {}".format(number_of_algorithms))

        if code or derivation or experiment_code:
            directory_name = os.path.join(config.output_code_path, output_name)
            if not os.path.exists(directory_name):
                os.makedirs(directory_name)
            else:
                # Removing existing algorithm files
                algorithms_dir_name = os.path.join(directory_name, config.language.name, subdir_name)
                cgu.remove_files(algorithms_dir_name)
                derivation_dir_name = os.path.join(directory_name, config.language.name, subdir_name, "derivation")
                cgu.remove_files(derivation_dir_name)

        if code or derivation:
            for n, (path, cost) in enumerate(paths):
            
                kernels_and_equations = []
                current_node = self.root
                for idx in path:
                    edge_label = current_node.edge_labels[idx]
                    kernels_and_equations.append(current_node.original_equations[idx])
                    kernels_and_equations.extend(edge_label.matched_kernels)
                    current_node = current_node.successors[idx]

                algorithm = cgu.Algorithm(self.input, current_node.equations, kernels_and_equations, cost)

                if code:
                    cgu.algorithm_to_file(output_name, subdir_name, algorithm_name.format(n), algorithm.code(), algorithm.experiment_input, algorithm.experiment_output)

                if derivation:
                    file_name = os.path.join(config.output_code_path, output_name, config.language.name, subdir_name, "derivation", "algorithm{}.txt".format(n))
                    directory_name = os.path.dirname(file_name)
                    if not os.path.exists(directory_name):
                        os.makedirs(directory_name)
                    output_file = open(file_name, "wt")
                    output_file.write(algorithm.derivation())
                    output_file.close()
                    if config.verbosity >= 2:
                        print("Generate derivation file {}".format(file_name))

        if experiment_code:

            reference_code.generate_reference_code(output_name, self.input)

            operand_generation.generate_operand_generator(output_name, self.input)

            algorithms = [(subdir_name, algorithm_name.format(0))]
            runner.generate_runner(output_name, algorithms)

        return number_of_algorithms

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


    def metric_mins(self):
        """Computes the minimal metric for each metric class.

        Returns a dictionary that maps the class (i.e. the first value in the
        metric) to the minimum of the second value of that class.
        """
        mins = dict()
        for node in self.active_nodes:
            if node.metric is None:
                node.metric = node.equations.metric()
            try:
                min, cost = mins[node.metric[0]]
            except KeyError:
                mins[node.metric[0]] = (node.metric[1], node.accumulated_cost)
            else:
                if min > node.metric[1]:
                    mins[node.metric[0]] = (node.metric[1], node.accumulated_cost)

        return mins

    def DS_tricks(self):

        new_nodes = []

        for node in self.active_nodes:
            transformed = []

            for equations in generate_variants(node.equations):
                transformed.extend(self.TR_tricks(equations))

            new_nodes.extend(self.create_nodes(node, *transformed))

        return new_nodes

    def DS_CSE_replacement(self):
        """Replaces common subexpression.

        Creates a new active node for each CSE. Current active nodes are not
        removed, such that the original equations are used for futher
        derivations as well.
        """

        new_nodes = []

        for node in self.active_nodes:
            if DS_step.CSE in node.applied_DS_steps:
                continue

            transformed = []

            for equations in generate_variants(node.equations):
                transformed.extend(self.TR_CSE_replacement(equations))

            new_nodes.extend(self.create_nodes(node, *transformed, previous_DS_step=DS_step.CSE))

        return new_nodes

    def DS_prune(self, mins):
        if self.active_nodes == [self.root]:
           return 0

        pruned_nodes = []
        for node in self.active_nodes:
            min, cost = mins[node.metric[0]]
            if node.metric[1] > (min + 4) and not node.successors:
                pruned_nodes.append(node)

        for node in pruned_nodes:
            node.labels.append("pruned")
            self.active_nodes.remove(node)

        return len(pruned_nodes)


    def TR_tricks(self, equations):

        for eqn_idx, equation in enumerate(equations):
            if not isinstance(equation.rhs, Symbol):
                for expr, position in equation.rhs.preorder_iter():
                    for callback_func, substitution in tricks.trick_MA.match(expr):
                        yield callback_func(substitution, equations, eqn_idx, position)
                # This break ensures that tricks are only applied to the first non-trivial equation.
                # Without it, unnecessary diamonds can appear in the derivation graph.
                break


    def TR_matrix_chain(self, equations, eqn_idx, initial_pos, explicit_inversion=False):

        try:
            msc = mcs.MatrixChainSolver(equations[eqn_idx][initial_pos], explicit_inversion)
        except mcs.MatrixChainNotComputable:
            return

        replacement = msc.tmp
        matched_kernels = msc.matched_kernels

        new_equation = matchpy.replace(equations[eqn_idx], initial_pos, replacement)
        equations_copy = equations.set(eqn_idx, new_equation)
        equations_copy = equations_copy.to_normalform()

        temporaries.set_equivalent_upwards(equations[eqn_idx].rhs, equations_copy[eqn_idx].rhs)
        
        yield (equations_copy, matched_kernels, equations)


    def TR_addition(self, equations, eqn_idx, initial_pos):

        # Note: For addition, we decided to only use a binary kernel, no
        #       variadic addition.

        new_expr, matched_kernels = decompose_sum(equations[eqn_idx][initial_pos])

        new_equation = matchpy.replace(equations[eqn_idx], initial_pos, new_expr)
        equations_copy = equations.set(eqn_idx, new_equation)
        equations_copy = equations_copy.to_normalform()

        yield (equations_copy, matched_kernels, equations)


    def TR_CSE_replacement(self, equations):

        for new_equations in CSEs.find_CSEs(equations):
            yield (new_equations, (), equations)


    def TR_unary_kernels(self, equations, eqn_idx, initial_pos):

        initial_node = equations[eqn_idx][initial_pos]

        # iterate over all subexpressions
        for node, _pos in initial_node.preorder_iter():
            pos = initial_pos + _pos

            kernel, substitution = select_optimal_match(collections_module.unary_kernel_DN.match(node))

            if kernel:
                # replacement
                
                matched_kernel = kernel.set_match(substitution, False)
                evaled_repl = matched_kernel.replacement

                # replace node with modified expression

                equations_copy = equations.set(eqn_idx, matchpy.replace(equations[eqn_idx], pos, evaled_repl))
                equations_copy = equations_copy.to_normalform()

                yield (equations_copy, (matched_kernel,), equations)


    def DS_factorizations(self):
        """Factorization derivation step.

        This derivation step applies factorizations to the active nodes. For a
        description of how exactly factorizations are applied, see the docstring
        of TR_factorizations().
        """

        new_nodes = []

        for node in self.active_nodes:
            if DS_step.factorizations in node.applied_DS_steps:
                continue

            # find all matrices that need to factored
            operands_to_factor = find_operands_to_factor(node.equations)

            operands_to_factor -= node.factored_operands
            if operands_to_factor:
                # construct dict {operand name: list of matched kernels for all valid factorizations}
                # It is constructed here because it can be reused. 
                factorization_dict = dict()
                for operand in operands_to_factor:
                    for type in collections_module.factorizations_by_type:
                        for kernel in type:
                            matches = list(matchpy.match(operand, kernel.pattern))

                            # This is important. Otherwise, the "break" at the end of
                            # the loop is reached. In that case, only the first
                            # factorization of each type is applied.
                            if not matches:
                                continue

                            # this is kind of stupid, there is only one match
                            for match_dict in matches:
                                matched_kernel = kernel.set_match(match_dict, False)
                                factorization_dict.setdefault(operand.name, []).append(matched_kernel)
                            break

                transformed = []

                found_symbol_inv_occurrence = []
                for equations in generate_variants(node.equations):
                    transformed.extend(self.TR_factorizations(equations, operands_to_factor, factorization_dict))

                # node.add_factored_operands(operands_to_factor)
                new_nodes.extend(self.create_nodes(node, *transformed, factored_operands=operands_to_factor, previous_DS_step=DS_step.factorizations))

        return new_nodes


    def TR_factorizations(self, equations, operands_to_factor, factorization_dict):
        """This function generates new equations by applying factorizations.

        For applying factorizations, a number of rules are applied:
        - Matrices are only factored if they have the ADMITS_FACTORIZATION
          property (this is tested by find_operands_to_factor() ).
        - We never apply different facotrizations to differrent occurrences of
          the same matrix. As an example, if the matrix A shows up twice, we
          will never apply LU to one occurrence and QR to another.
        - Factorization will always be applied to matrices that appear
          immediately inside an inverse. That is, the A in Inverse(A) will be
          factored. A and B in Inverse(Times(A, B)) don't have to be factored.
        - If there is a summand that contains a matrix which is factored, but
          the summand does not contain any occurrences of that matrix within an
          inverse, that matrix will not be factored in this summand. This is
          to make sure that for A+inv(A), the first A is not factored.
        - Some factorizations rule out others. If Cholesky can be applied, LU
          will no be applied. Which factorization can be applied per operand is
          decided in DS_factorizations(). The factorization_dict contains the
          valid factorizations.

        This function applies factorization in all possible combinations that
        obey the rules above.
        """

        transformed_expressions = []

        # find all occurrences
        all_occurrences = list(find_occurrences(equations, operands_to_factor))

        blocking_products = list(find_blocking_products(equations, operands_to_factor))

        # Removing groups (summands) which do not contain any inverted occurrences.
        candidate_occurrences = []
        for oc_group in group_occurrences(all_occurrences):
            if any(oc.type != InverseType.none for oc in oc_group):
                candidate_occurrences.extend(oc_group)

        # collect all operands that show up
        ops = set(oc.operand.name for oc in candidate_occurrences)

        # Symbols directely inside an inverse always have to be factored.
        ops_must_factor = set()
        ops_may_factor = set()
        for op in ops:
            if any(oc.operand.name == op and oc.symbol for oc in candidate_occurrences):
                ops_must_factor.add(op)
            else:
                ops_may_factor.add(op)

        for ops_subset in powerset(ops_may_factor):

            factor_ops = ops_must_factor.union(ops_subset)
            if not factor_ops or factor_ops in blocking_products:
                continue

            factorizations_candidates = []
            for op in factor_ops:
                factorizations_candidates.append(factorization_dict[op])

            # apply all factorizations
            for factorizations in itertools.product(*factorizations_candidates):
                facts_dict = dict(zip(factor_ops, factorizations))

                # collect matched kernels (avoiding duplicates)
                matched_kernels = []
                _already_seen = set()
                for matched_kernel in factorizations:
                    if matched_kernel.id not in _already_seen:
                        matched_kernels.append(matched_kernel)
                        _already_seen.add(matched_kernel.id)

                # collect replacements 
                replacements_per_equation = dict()

                for oc in candidate_occurrences:
                    if oc.operand.name in factor_ops:
                        replacements_per_equation.setdefault(oc.eqn_idx, []).append((oc.position, facts_dict[oc.operand.name].replacement))

                # replace
                equations_list = list(equations.equations)

                for eqn_idx, replacements in replacements_per_equation.items():
                    if replacements:
                        equations_list[eqn_idx] = matchpy.replace_many(equations_list[eqn_idx], replacements)
                
                equations_copy = Equations(*equations_list)
                equations_copy = equations_copy.simplify()
                equations_copy.set_equivalent(equations)
                equations_copy = equations_copy.to_SOP().simplify()

                transformed_expressions.append((equations_copy, matched_kernels, equations))


        return transformed_expressions


    def init_temporaries(self, equations):

        seen_before = set()
        operands_to_factor = find_operands_to_factor(equations)
        for equations_var in generate_variants(equations):
            for equation in equations_var:
                for inv_expr, pos in find_explicit_symbol_inverse(equation.rhs):
                    if inv_expr not in seen_before:
                        special_properties.add_expression(inv_expr, [])
                        seen_before.add(inv_expr)
                for expr, pos in equation.rhs.preorder_iter():
                    # if isinstance(expr, Times) and any(is_inverse(operand) and (not isinstance(operand.operand, Symbol) or operand.operand in operands_to_factor) for operand in expr.operands):
                    if (isinstance(expr, Times) and any(is_inverse(operand) and operand.operand in operands_to_factor for operand in expr.operands)) or (is_inverse(expr) and not isinstance(expr.operand, Symbol)):
                        # TODO is it dangerous to add inv(expr) here? It becomes explicit inversion, even if it originally wasn't.
                        if expr not in seen_before:
                            special_properties.add_expression(expr, [])
                            seen_before.add(expr)

            # TODO check_consistency for Equations?
            check_consistency(equation)

class DerivationGraphNode(base.GraphNodeBase):

    _counter = 0

    def __init__(self, equations=None, predecessor=None, factored_operands=None, previous_DS_step=None):
        super().__init__(predecessor, factored_operands, previous_DS_step)

        # IDs for dot output
        self.id = DerivationGraphNode._counter
        self.name = "".join(["node", str(self.id)])
        DerivationGraphNode._counter +=1
        self.equations = equations

        
    def get_payload(self):
        return self.equations

    def is_terminal(self):
        """Returns true if the expression of self is fully computed."""
        return all([isinstance(equation.rhs, Symbol) for equation in self.equations])
