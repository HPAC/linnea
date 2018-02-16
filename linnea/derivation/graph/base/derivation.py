from . import base

from ..utils import generate_variants, find_operands_to_factor, \
                    find_occurrences, InverseType, powerset, group_occurrences

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
from ...matrix_sum import decompose_sum

from ...utils import select_optimal_match
from ....code_generation.memory import memory as memory_module
from ....code_generation import utils as cgu
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


    def create_nodes(self, predecessor, *description, factored_operands=None):
        new_nodes = []
        # print(description)
        # if description:

        _factored_operands = predecessor.factored_operands
        if factored_operands:
            _factored_operands = _factored_operands.union(factored_operands)
            
        for equations, metric, edge_label in description:
            new_node = DerivationGraphNode(equations, metric, predecessor, _factored_operands)
            self.nodes.append(new_node)
            new_nodes.append(new_node)
            predecessor.set_labeled_edge(new_node, edge_label)
        return new_nodes


    def write_output(self, code=True, pseudocode=False, output_name="tmp", operand_generator=False, algorithms_limit=1, graph=False):

        if not config.output_path:
            raise config.OutputPathNotSet("Unable to write output: output_path not set.")

        if graph:
            self.write_graph(output_name)
        
        algorithm_paths = list(self.all_algorithms())
        algorithm_paths.sort(key=operator.itemgetter(1))

        number_of_algorithms = len(algorithm_paths)
        if config.verbosity >= 1:
            self.print("Number of algorithms: {}".format(number_of_algorithms))
        if number_of_algorithms > algorithms_limit:
           algorithm_paths = algorithm_paths[:algorithms_limit]

        #if number_of_algorithms == 0:
        #    print("No algorithm generated for this example")
        #    return False

        if code or pseudocode or operand_generator:
            directory_name = os.path.join(config.output_path, output_name)
            if not os.path.exists(directory_name):
                os.makedirs(directory_name)

        # TODO what is the purpose of "len(algorithm_paths[0][0]) > 0"?
        # if number_of_algorithms > 0 and len(algorithm_paths[0][0]) > 0:
        for n, (algorithm_path, cost) in enumerate(algorithm_paths):
        
            matched_kernels = []
            current_node = self.root
            for idx in algorithm_path:
                edge_label = current_node.edge_labels[idx]
                matched_kernels.extend(edge_label.matched_kernels)
                current_node = current_node.successors[idx]

            algorithm = cgu.Algorithm(self.root.equations, current_node.equations, matched_kernels, cost)

            if code:
                # TODO change
                # here, I need to generat functions
                cgu.algorithm_to_file(output_name, "algorithm{}".format(n), algorithm.code(), algorithm.experiment_input, algorithm.experiment_output)
                # code_gen.algorithm_to_file(output_name, algorithm_name, algorithm, input, output)
                # file_name = os.path.join(config.output_path, config.language.name, output_name, "algorithms", "algorithm{}{}".format(str(n), config.filename_extension))
                # output_file = open(file_name, "wt")
                # # try:
                # #     output_file.write(algorithm.code())
                # # except KeyError:
                # #     pass
                # _code = algorithm.code()
                # output_file.write(_code)
                # output_file.close()

            if pseudocode:
                file_name = os.path.join(config.output_path, output_name, config.language.name, "pseudocode", "algorithm{}.txt".format(n))
                directory_name = os.path.dirname(file_name)
                if not os.path.exists(directory_name):
                    os.makedirs(directory_name)
                output_file = open(file_name, "wt")
                output_file.write(algorithm.pseudocode())
                output_file.close()

        if operand_generator:
            input, output = self.root.equations.input_output()
            input_str = ", ".join([operand.name for operand in input])
            output_str = ", ".join([operand.name for operand in output])
            cgu.algorithm_to_file(output_name, "naive", self.root.equations.to_julia_expression(), input_str, output_str, config.Language.Julia)
            cgu.algorithm_to_file(output_name, "recommended", self.root.equations.to_julia_expression(recommended=True),
                                  input_str, output_str, config.Language.Julia)
            cgu.algorithm_to_file(output_name, "naive", self.root.equations.to_cpp_expression(config.CppLibrary.Blaze),
                                  input_str, output_str, config.Language.Cpp, ".hpp", "blaze")
            cgu.algorithm_to_file(output_name, "naive", self.root.equations.to_cpp_expression(config.CppLibrary.Eigen),
                                  input_str, output_str, config.Language.Cpp, ".hpp", "eigen")
            cgu.algorithm_to_file(output_name, "naive", self.root.equations.to_cpp_expression(config.CppLibrary.Armadillo),
                                  input_str, output_str, config.Language.Cpp, ".hpp", "armadillo")
            cgu.algorithm_to_file(output_name, "recommended",
                                  self.root.equations.to_cpp_expression(config.CppLibrary.Eigen, recommended=True),
                                  input_str, output_str, config.Language.Cpp, ".hpp", "eigen")
            cgu.algorithm_to_file(output_name, "recommended",
                                  self.root.equations.to_cpp_expression(config.CppLibrary.Armadillo, recommended=True),
                                  input_str, output_str, config.Language.Cpp, ".hpp", "armadillo")
            cgu.algorithm_to_file(output_name, "naive", self.root.equations.to_julia_expression(), input_str, output_str,
                                  config.Language.Matlab, ".m")
            cgu.algorithm_to_file(output_name, "recommended", self.root.equations.to_julia_expression(recommended=True),
                                  input_str, output_str, config.Language.Matlab, ".m")
            cge.operand_generator_to_file(output_name, input, input_str)
            cge.operand_generator_to_file(output_name, input, input_str, language=config.Language.Cpp)
            cge.operand_generator_to_file(output_name, input, input_str, language=config.Language.Matlab)
            cge.benchmarker_to_file(output_name, algorithms_count=len(algorithm_paths), language=config.Language.Julia)
            cge.benchmarker_to_file(output_name, language=config.Language.Matlab)
            cge.benchmarker_to_file(output_name, language=config.Language.Cpp)

            # create language runner

        # TODO missing
        # - change paths and use name
        return True

    def optimal_algorithm_to_str(self):
        matched_kernels, cost, final_equations = self.optimal_algorithm()
        if matched_kernels:
            algorithm = cgu.Algorithm(self.root.equations, final_equations, matched_kernels, cost)
            code = algorithm.code()
            # code_lines = code.splitlines()
            # del code_lines[1:3] # remove "using Base.LinAlg..."
            # code = "\n".join(code_lines)
            function_name = "linnea_function{:X}".format(hash(self.root.equations))
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


    def TR_matrix_chain(self, equations, eqn_idx, initial_pos, metric, explicit_inversion=False):
        equations_copy = copy.deepcopy(equations)
        expr = equations_copy[eqn_idx][initial_pos]

        try:
            #print("before")
            msc = mcs.MatrixChainSolver(expr, explicit_inversion)
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


    def DS_factorizations(self):

        new_nodes = []
        inactive_nodes = []

        # store list of matrices that have been factored already

        # just_factored_operands = set()

        for node in self.active_nodes:


            # find all matrices that need to factored
            operands_to_factor = find_operands_to_factor(node.equations)
            # print(operands_to_factor)

            # just_factored_operands.update(operands_to_factor)
            operands_to_factor -= node.factored_operands
            # print(operands_to_factor)
            if operands_to_factor:
                # construct dict {operand name: list of matched kernels for all valid factorizations}
                factorization_dict = dict()
                for operand in operands_to_factor:
                    for type in collections_module.factorizations_by_type:
                        for kernel in type:
                            # print()
                            matches = list(matchpy.match(operand, kernel.pattern))

                            # This is important. Otherwise, the "break" at the end of
                            # the loop is reached. In that case, only the first
                            # factorization of each type is applied.
                            if not matches:
                                continue

                            # this is kind of stupid, there is only one match
                            for match_dict in matches:
                                matched_kernel = kernel.set_match(match_dict, False, CSE_rules=True)
                                factorization_dict.setdefault(operand.name, []).append(matched_kernel)
                            break
                # print(factorization_dict)

                transformed = []

                found_symbol_inv_occurrence = False
                for equations in generate_variants(node.equations):
                    _transformed, _found_symbol_inv_occurrence = self.TR_factorizations(equations, operands_to_factor, factorization_dict)
                    found_symbol_inv_occurrence = found_symbol_inv_occurrence or _found_symbol_inv_occurrence
                    transformed.extend(_transformed)

                new_nodes.extend(self.create_nodes(node, *transformed, factored_operands=node.factored_operands.union(operands_to_factor)))

                # if there are no symbol inverses, the current node stays active
                # because it's possible to make progress without factorizations
                if found_symbol_inv_occurrence:
                    inactive_nodes.append(node)
        
        for node in inactive_nodes:
            self.active_nodes.remove(node)

        self.add_active_nodes(new_nodes)

        return len(new_nodes)



    def TR_factorizations(self, equations, operands_to_factor, factorization_dict):

        transformed_expressions = []

        # find all occurrences
        all_occurrences = list(find_occurrences(equations, operands_to_factor))

        # print(all_occurrences)

        # filter out occurrences that have InverseType.none
        # filter our occurrences with symbol=True
        non_inv_occurrences = []
        inv_occurrences = []
        symbol_inv_occurrences = []
        found_symbol_occurrence = False
        for occurrence in all_occurrences:
            if occurrence.symbol:
                symbol_inv_occurrences.append(occurrence)
                found_symbol_occurrence = True
            elif occurrence.type is InverseType.none:
                non_inv_occurrences.append(occurrence)
            else:
                inv_occurrences.append(occurrence)

        # print(non_inv_occurrences)

        # group InverseType.none occurrences by operands
        non_inv_occurrences_dict = dict()
        for occurrence in non_inv_occurrences:
            non_inv_occurrences_dict.setdefault(occurrence.operand.name, []).append(occurrence)

        # print(non_inv_occurrences_dict)

        inv_occurrences_grouped = group_occurrences(inv_occurrences)
        symbol_inv_occurrences_grouped = group_occurrences(symbol_inv_occurrences)

        # collect all operands that show up
        ops = set()

        for ocs in symbol_inv_occurrences_grouped:
            ops.update([oc.operand.name for oc in ocs])

        # print(inv_occurrences_grouped)
        # all subsets of grouped occurrences
        for subset_inv_occurrences_grouped in powerset(inv_occurrences_grouped):

            # print(subset_inv_occurrences_grouped)

            ops_copy = ops.copy()
            for ocs in subset_inv_occurrences_grouped:
                ops_copy.update([oc.operand.name for oc in ocs])

            ops_copy = list(ops_copy)
            # print("OPS", ops_copy)
            
            # collect all non inverse occurences for those operands
            non_inv_occurrences = []
            factorizations_candidates = []
            for op in ops_copy:
                try:
                    ocs = non_inv_occurrences_dict[op]
                except KeyError:
                    pass
                else:
                    non_inv_occurrences.extend(ocs)
                
                factorizations_candidates.append(factorization_dict[op])

            # print(factorizations_candidates)
            # print(non_inv_occurrences)
            if not factorizations_candidates:
                continue

            # generated all subsets
            for subset_non_inv_occurrences in powerset(non_inv_occurrences):
                # ocs = subset_inv_occurrences_grouped + subset_non_inv_occurrences

                # TODO do we allow different factorizations_candidates for the same operand?
                # factorizations_candidates = []
                # for oc in subset_inv_occurrences_grouped:
                #     factorizations_candidates.append(factorization_dict[oc[0][2].name])

                # for oc in subset_non_inv_occurrences:
                #     factorizations_candidates.append(factorization_dict[oc[2].name])


                # apply all factorizations
                for factorizations in itertools.product(*factorizations_candidates):
                    facts_dict = dict(zip(ops_copy, factorizations))

                    # print(facts_dict)

                    # TODO missing: label explicit inversions as such

                    # collect matched kernels (avoiding duplicates)
                    matched_kernels = []
                    _already_seen = set()
                    for matched_kernel in factorizations:
                        if matched_kernel.id not in _already_seen:
                            matched_kernels.append(matched_kernel)
                            _already_seen.add(matched_kernel.id)

                    # collect replacements 
                    replacements_per_equation = dict()

                    for ocs in symbol_inv_occurrences_grouped:
                        for oc in ocs:
                            replacements_per_equation.setdefault(oc.eqn_idx, []).append((oc.position, facts_dict[oc.operand.name].replacement))

                    for ocs in subset_inv_occurrences_grouped:
                        for oc in ocs:
                            replacements_per_equation.setdefault(oc.eqn_idx, []).append((oc.position, facts_dict[oc.operand.name].replacement))

                    for oc in subset_non_inv_occurrences:
                        replacements_per_equation.setdefault(oc.eqn_idx, []).append((oc.position, facts_dict[oc.operand.name].replacement))

                    # replace
                    equations_copy = copy.deepcopy(equations)

                    for eqn_idx, replacements in replacements_per_equation.items():
                        if replacements:
                            equations_copy[eqn_idx] = matchpy.replace_many(equations_copy[eqn_idx], replacements)
                    
                    equations_copy.to_normalform()
                    equations_copy.set_equivalent(equations)

                    # print(equations_copy)
                    # print(matched_kernels)                   

                    edge_label = base.EdgeLabel(*matched_kernels)
                    transformed_expressions.append((equations_copy, equations_copy.metric(), edge_label))
                        
        return transformed_expressions, found_symbol_occurrence


class DerivationGraphNode(base.GraphNodeBase):

    _counter = 0

    def __init__(self, equations=None, metric=None, predecessor=None, factored_operands=None):
        super(DerivationGraphNode, self).__init__(metric, predecessor)

        # IDs for dot output
        self.id = DerivationGraphNode._counter
        self.name = "".join(["node", str(self.id)])
        DerivationGraphNode._counter +=1

        self.equations = equations
        if factored_operands is None:
            self.factored_operands = set()
        else:
            self.factored_operands = factored_operands
        
    def get_payload(self):
        return self.equations

    def is_terminal(self):
        """Returns true if the expression of self is fully computed."""
        return all([isinstance(equation.rhs, Symbol) for equation in self.equations])
