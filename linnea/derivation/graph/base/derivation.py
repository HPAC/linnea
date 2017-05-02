from . import base

from ..utils import generate_variants

from ....algebra.expression import Symbol, ConstantScalar, \
                                   Operator, Plus, Times, Equal, Transpose

from ....algebra.transformations import simplify
from ....algebra.representations import to_SOP
from ....algebra.properties import Property as properties

import matchpy
import itertools
import copy
import operator
import os.path

from ... import tricks
from ... import CSEs
from ... import matrix_chain_solver as mcs
# from ... import matrix_chain_solver_new as mcs
from ...matrix_sum import decompose_sum

from ...utils import select_optimal_match
from ....code_generation.memory import memory as memory_module
from ....code_generation.utils import Algorithm
from ....code_generation import experiments as cge
from .... import config
from .... import temporaries

collections_module = config.import_collections()

class DerivationGraphBase(base.GraphBase):
    """

    Graph for equations.
    """

    def __init__(self, root_equations):
        super(DerivationGraphBase, self).__init__()
        self.root = DerivationGraphNode(root_equations)
        self.active_nodes = [self.root]
        self.nodes = [self.root]
        self.level_counter = -1


    def create_nodes(self, predecessor, *description):
        new_nodes = []
        # print(description)
        # if description:
        for equations, metric, edge_label in description:
            new_node = DerivationGraphNode(equations, metric, predecessor)
            self.nodes.append(new_node)
            new_nodes.append(new_node)
            predecessor.set_labeled_edge(new_node, edge_label)
        return new_nodes


    def algorithms_to_files(self, pseudocode=False, experiments=False, name=None):
        # TODO additional optional argument all/optimal to specify whether to
        # generate all algorithms or only the optimal one

        algorithm_paths = list(self.all_algorithms(self.root))
        algorithm_paths.sort(key=operator.itemgetter(1))

        # algorithm_paths = [self.optimal_algorithm_path()]

        number_of_algorithms = len(algorithm_paths)
        self.print("Number of algorithms: {}".format(number_of_algorithms))
        if number_of_algorithms > 100:
            algorithm_paths = algorithm_paths[0:100]

        for n, (algorithm_path, cost) in enumerate(algorithm_paths):
            
            matched_kernels = []
            current_node = self.root
            for idx in algorithm_path:
                edge_label = current_node.edge_labels[idx]
                matched_kernels.extend(edge_label.matched_kernels)
                current_node = current_node.successors[idx]

            algorithm = Algorithm(self.root.equations, current_node.equations, matched_kernels, cost)
            # algorithm.liveness_analysis()
            # algorithm.pseudocode()

            file_name = os.path.join(config.output_path, config.language.name, "algorithms", "algorithm{}{}".format(str(n), config.filename_extension))
            output_file = open(file_name, "wt")
            # try:
            #     output_file.write(algorithm.code())
            # except KeyError:
            #     pass
            code = algorithm.code()
            output_file.write(code)
            output_file.close()

            if pseudocode:
                file_name = os.path.join(config.output_path, "pseudocode", "algorithm{}.txt".format(n))
                output_file = open(file_name, "wt")
                output_file.write(algorithm.pseudocode())
                output_file.close()

            if experiments:
                code_lines = code.splitlines()
                del code_lines[1:3] # remove "using Base.LinAlg..."
                code = "\n".join(code_lines)
                cge.experiment_to_file(name, "algorithm{}".format(n), code, algorithm.experiment_input, algorithm.experiment_output)

        if experiments:
            input, output = self.root.equations.input_output()
            input_str = ", ".join([operand.name for operand in input])
            output_str = ", ".join([operand.name for operand in output])
            cge.experiment_to_file(name, "naive", self.root.equations.to_julia_expression(), input_str, output_str)
            cge.operand_generator_to_file(name, input, input_str)


    def metric_mins(self):
        """Computes the minimal metric for each metric class.

        Returns a dictionary that maps the class (i.e. the first value in the
        metric) to the minimum of the second value of that class.
        """
        mins = dict()
        for node in self.active_nodes:
            try:
                min = mins[node.metric[0]]
            except KeyError:
                mins[node.metric[0]] = node.metric[1]
            else:
                if min > node.metric[1]:
                    mins[node.metric[0]] = node.metric[1]

        return mins

    def DS_tricks(self):

        new_nodes = []

        for node in self.active_nodes:
            transformed = []

            for equations in generate_variants(node.equations):
                # new_nodes.extend(self.create_nodes(node, self.TR_tricks(equations)))
                transformed.extend(self.TR_tricks(equations))

            new_nodes.extend(self.create_nodes(node, *transformed))

        self.add_active_nodes(new_nodes)
        # self.active_nodes.extend(new_nodes)

        return len(new_nodes)

    def DS_CSE_replacement(self):
        """Replaces common subexpression.

        Creates a new active node for each CSE. Current active nodes are not
        removed, such that the original equations are used for futher
        derivations as well.
        """

        ###############
        # Replace CSEs.

        new_nodes = []

        for node in self.active_nodes:

            transformed = []

            for equations in generate_variants(node.equations):
                transformed.extend(self.TR_CSE_replacement(equations))

            new_nodes.extend(self.create_nodes(node, *transformed))

        # We don't remove any active nodes here because for the derivation, we
        # also want to continue with the original equations were no CSEs where
        # replaced.

        self.add_active_nodes(new_nodes)
        # self.active_nodes.extend(new_nodes)

        return len(new_nodes)

    def DS_prune(self, mins):
        if self.active_nodes == [self.root]:
           return 0

        pruned_nodes = []
        for node in self.active_nodes:
            min = mins[node.metric[0]]
            if node.metric[1] > (min + 4):
                pruned_nodes.append(node)
                for predecessor in node.predecessors:
                    predecessor.remove_edge(node)

        self.remove_nodes(pruned_nodes)

        return len(pruned_nodes)



    def TR_tricks(self, equations):

        transformed_expressions = []

        for eqn_idx, equation in enumerate(equations):
            if not isinstance(equation.rhs, Symbol):
                for node, position in equation.rhs.preorder_iter():
                    for callback_func, substitution in tricks.trick_MA.match(node):
                        transformed_expressions.append(callback_func(substitution, equations, eqn_idx, position))
                # Ensures that tricks are only applied to the first non-trivial
                # equation. Without this break, unnecessary diamonds can appear
                # in the derivation graph.
                break

        return transformed_expressions


    def TR_matrix_chain(self, equations, eqn_idx, initial_pos, metric):
        equations_copy = copy.deepcopy(equations)
        expr = equations_copy[eqn_idx][initial_pos]

        try:
            # print("before")
            msc = mcs.MatrixChainSolver(expr)
        except mcs.MatrixChainNotComputable:
            return []

        replacement = msc.tmp
        matched_kernels = msc.matched_kernels

        # print("#####")
        # print(equations_copy[eqn_idx].get_successor(initial_pos[:-1]))
        # print(equations)
        # print(expr)
        # print(matched_kernels)
        # print(replacement)
        # print(eqn_idx, initial_pos)
        
        equations_copy[eqn_idx] = matchpy.replace(equations_copy[eqn_idx], initial_pos, replacement)

        equations_copy[eqn_idx] = to_SOP(simplify(equations_copy[eqn_idx]))

        # # equivalent expression experiment
        # path_before, path_after = temporaries.expr_diff(equations[eqn_idx].rhs, equations_copy[eqn_idx].rhs)
        # # print(path_before, path_after)
        # expr_before = equations[eqn_idx].rhs[path_before]
        # expr_after = equations_copy[eqn_idx].rhs[path_after]
        # tmp = temporaries.create_tmp(expr_before, True)
        # # print(tmp)
        # print("added in MC: ", expr_before, temporaries._get_equivalent(expr_after), tmp)
        # temporaries._table_of_temporaries[str(temporaries._get_equivalent(expr_after))] = tmp

        temporaries.set_equivalent_upwards(equations[eqn_idx].rhs, equations_copy[eqn_idx].rhs)

        new_metric = equations_copy.metric()
        if new_metric <= metric:
            # print(equations_copy)
            return [(equations_copy, new_metric, base.EdgeLabel(*matched_kernels))]
        else:
            return []


    def TR_addition(self, equations, eqn_idx, initial_pos, metric):
        equations_copy = copy.deepcopy(equations)
        expr = equations_copy[eqn_idx][initial_pos]

        # Note: For addition, we decided to only use a binary kernel, no
        #       variadic addition.

        expr, matched_kernels = decompose_sum(expr)

        equations_copy[eqn_idx] = matchpy.replace(equations_copy[eqn_idx], initial_pos, expr)

        equations_copy[eqn_idx] = to_SOP(simplify(equations_copy[eqn_idx]))

        new_metric = equations_copy.metric()
        if new_metric <= metric:
            return [(equations_copy, new_metric, base.EdgeLabel(*matched_kernels))]
        else:
            return []


    def TR_CSE_replacement(self, equations):

        # for plus
        sums = []
        eqn_indices_plus = []
        sum_positions = []

        # for times
        products = []
        eqn_indices_times = []
        product_positions = []

        # for general
        expressions = []
        expr_positions = []

        for n, equation in enumerate(equations):
            for node, pos in equation.preorder_iter():
                # Plus
                if isinstance(node, Plus):
                    sums.append(node)
                    eqn_indices_plus.append(n)
                    sum_positions.append((n, pos))
                # Times
                if isinstance(node, Times):# and CSEs._is_simple_times_CSE(node):
                    products.append(node)
                    eqn_indices_times.append(n)
                    product_positions.append((n, pos))
                # General
                if not isinstance(node, Times) and not isinstance(node, Symbol) and not isinstance(node, ConstantScalar) and not isinstance(node, Equal) and not (isinstance(node, Operator) and node.arity is matchpy.Arity.unary and isinstance(node.operand, Symbol)):
                    # Products, symbols and Unary(Symbol) are not considered.
                    expressions.append(node)
                    expr_positions.append((n, pos))

        transformed_expressions = []
        if sums:
            transformed_expressions.extend(CSEs.CSE_replacement_plus(equations, sums, sum_positions))
        if products:
            transformed_expressions.extend(CSEs.CSE_replacement_times(equations, products, product_positions))
        if expressions:
            transformed_expressions.extend(CSEs.CSE_replacement_general(equations, expressions, expr_positions))

        return transformed_expressions


    # @profile
    def TR_factorizations(self, equations, eqn_idx, initial_pos, metric):
        """

        This function takes priorities into
        account. That is, for a given matrix, from each 1, 2 and 3, the best
        factorization is chosen (the best one is the first one that is
        applicable).
        1. Cholesky, LDL, PLU
        2. QR
        3. Eigen, SVD
        """
        transformed_expressions = []        

        initial_node = equations[eqn_idx][initial_pos]

        known_node = []
        for node, _pos in initial_node.preorder_iter():
            full_pos = copy.copy(initial_pos)
            full_pos.extend(_pos)
            if node not in known_node:
                known_node.append(node)
            else:
                continue

            # TODO is this necessary?
            # if not node.has_property(properties.ADMITS_FACTORIZATION):
            #     continue

            for type in collections_module.factorizations_by_type:
                for kernel in type:
                    # for each match, generate copy and apply kernel
                    matches = list(matchpy.match(node, kernel.pattern))

                    # This is important. Otherwise, the "break" at the end of
                    # the loop is reached. In that case, only the first
                    # factorization of each type is applied.
                    if not matches:
                        continue

                    for match_dict in matches:
                        # print(kernel.signature, match_dict)
                        # replacement
                        equations_copy = copy.deepcopy(equations)

                        matched_kernel = kernel.set_match(match_dict, False, CSE_rules=True, blocked_products=True)
                        evaled_repl = matched_kernel.replacement

                        # replace node with modified expression
                        # print(initial_node[_pos[:-1]], evaled_repl)
                        equations_copy[eqn_idx] = matchpy.replace(equations_copy[eqn_idx], full_pos, evaled_repl)

                        # deal with additional occurrences of the replaced subexpression
                        common_subexp_rules = matched_kernel.CSE_rules

                        equations_copy.replace_all(common_subexp_rules)
                        
                        # # equivalent expression experiment
                        # # print("###")
                        # # print(equations[eqn_idx].rhs, equations_copy[eqn_idx].rhs)
                        # path_before, path_after = expr_diff(equations[eqn_idx].rhs, equations_copy[eqn_idx].rhs)
                        # # print(path_before, path_after)
                        # expr_before = equations[eqn_idx].rhs[path_before]
                        # expr_after = equations_copy[eqn_idx].rhs[path_after]
                        # print(temporaries._get_equivalent(expr_before), str(temporaries._get_equivalent(expr_before)) in temporaries._table_of_temporaries, temporaries._get_equivalent(expr_after))
                        # if str(temporaries._get_equivalent(expr_before)) in temporaries._table_of_temporaries:
                        #     print("stored: ", temporaries._table_of_temporaries[str(temporaries._get_equivalent(expr_before))])
                        # tmp = temporaries.create_tmp(expr_before, True)
                        # print(tmp)
                        # temporaries._table_of_temporaries[str(temporaries._get_equivalent(expr_after))] = tmp

                        temporaries.set_equivalent(equations[eqn_idx].rhs, equations_copy[eqn_idx].rhs)

                        new_metric = equations_copy.metric()
                        # print("metric", new_metric, metric)
                        # print(equations_copy)
                        if new_metric <= metric:
                            edge_label = base.EdgeLabel(matched_kernel)
                            transformed_expressions.append((equations_copy, new_metric, edge_label))

                    break
        return transformed_expressions


    def TR_unary_kernels(self, equations, eqn_idx, initial_pos, metric):

        transformed_expressions = []

        initial_node = equations[eqn_idx][initial_pos]

        # iterate over all subexpressions
        for node, _pos in initial_node.preorder_iter():
            pos = copy.copy(initial_pos)
            pos.extend(_pos)

            kernel, substitution = select_optimal_match(collections_module.unary_kernel_DN.match(node))

            if kernel:
                # replacement
                equations_copy = copy.deepcopy(equations)
                
                matched_kernel = kernel.set_match(substitution, False)
                evaled_repl = matched_kernel.replacement

                # replace node with modified expression

                equations_copy[eqn_idx] = matchpy.replace(equations_copy[eqn_idx], pos, evaled_repl)
                
                new_metric = equations_copy.metric()
                # print("metric", m, metric)
                # print(equations_copy)
                if new_metric <= metric:
                    edge_label = base.EdgeLabel(matched_kernel)
                    transformed_expressions.append((equations_copy, new_metric, edge_label))

        return transformed_expressions


class DerivationGraphNode(base.GraphNodeBase):

    _counter = 0

    def __init__(self, equations=None, metric=None, predecessor=None):
        super(DerivationGraphNode, self).__init__(metric, predecessor)

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
