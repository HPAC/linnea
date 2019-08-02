from .... import config

import itertools
import operator
import copy
import math
import os
import collections


class GraphBase():

    def __init__(self):
        super().__init__()
        self.step_counter = -1
        self.k_best_state = 1

    def print(self, str):
        if config.verbosity >= 1:
            print(str)

    def print_line(self):
        if config.verbosity >= 1:
            print("{:-<34}".format(""))

    def print_result_sep(self, str, val):
        if config.verbosity >= 1:
            print("{0:.<22}{1:.>12}".format(str, val))

    def print_result(self, str, val):
        if config.verbosity >= 1:
            print("{0:<22}{1:>12}".format(str, val))

    def derivation(self):
        raise NotImplementedError()

    def terminal_nodes(self):
        return list(filter(operator.methodcaller("is_terminal"), self.nodes))

    def nodes_count(self):
        return len(self.nodes)

    def create_nodes(self):
        raise NotImplementedError()

    def remove_node(self, node):
        self.nodes.remove(node)

    def remove_nodes(self, nodes):
        """removes nodes from the graph"""
        for node in nodes:
            self.remove_node(node)

    def to_dot(self, style):
        # TODO use source and sink ranks

        optimal_edges = set()

        if any(map(operator.methodcaller("is_terminal"), self.nodes)):

            optimal_path = next(self.k_shortest_paths(1))

            optimal_paths = list(itertools.takewhile(lambda path: path[1]==optimal_path[1], self.k_shortest_paths(100)))

            # print(minimal_cost)
            for path, cost in optimal_paths:
                current_node = self.root
                for path_idx in path:
                    previous_node_id = current_node.id
                    current_node = current_node.successors[path_idx]
                    optimal_edges.add((previous_node_id, current_node.id))

        out = ["digraph G {", "ranksep=2.5;", "rankdir=TB;"]
        out.extend([node.to_dot(style, optimal_edges) for node in self.nodes])
        out.append("}")
        return "\n".join(out)


    def write_graph(self, output_name, style, file_name="graph"):
        file_path = os.path.join(config.output_code_path, output_name, config.language.name, "{}.gv".format(file_name))
        directory_name = os.path.dirname(file_path)
        if not os.path.exists(directory_name):
            os.makedirs(directory_name)
        output_file = open(file_path, "wt")
        output_file.write(self.to_dot(style))
        output_file.close()
        if config.verbosity >= 2:
            print("Generate graph file {}".format(file_path))


    def all_algorithms(self):
        """Generates all paths in the graph.

        Careful: For large graphs, this is very slow. If not all paths are
            needed, use k_shortest_paths() instead.

        Yields:
            list: A path, in the form of indices of successors.
            float: The cost of the path.
        """
        stack = collections.deque([(self.root, [], 0)])
        while stack:
            node, path, cost = stack.popleft()
            if node.successors == [] and node.is_terminal():
                yield path, cost
            else:
                for idx, successor in enumerate(node.successors):
                    path_copy = copy.copy(path)
                    path_copy.append(idx)
                    stack.appendleft((successor, path_copy, cost + node.edge_labels[idx].cost))


    def shortest_path(self):
        try:
            path, cost = next(self.k_shortest_paths(1))
        except StopIteration:
            return ([], math.inf)
        else:
            return (path, cost)


    def optimal_algorithm(self):
        matched_kernels = []
        current_node = self.root
        path, cost = self.shortest_path()
        if path:
            for idx in path:

                matched_kernels.extend(current_node.edge_labels[idx].matched_kernels)
                if current_node.successors:
                    current_node = current_node.successors[idx]

            # beginning or end is missing
            return matched_kernels, cost, current_node.get_payload()
        else:
            return [], math.inf, None


    def topological_sort(self):
        stack = collections.deque()
        temp = set()
        perm = set()
        # TODO shouldn't it be sufficient to start at the root node?
        for node in self.nodes:
            node._topological_sort_visit(temp, perm, stack)
        return list(stack)


    def set_shortest_path_successor(self):
        sorted_nodes = self.topological_sort()

        cost_from_sink = dict()
        shortest_path_successor = dict() # TODO remove?
        for node in reversed(sorted_nodes):
            if not node.successors:
                if node.is_terminal():
                    cost_from_sink[node.id] = 0
                else:
                    cost_from_sink[node.id] = math.inf
            for predecessor in node.predecessors:
                idx = predecessor.successors.index(node)
                cost = predecessor.edge_labels[idx].cost

                if cost_from_sink.get(predecessor.id, math.inf) > cost_from_sink[node.id] + cost:
                    cost_from_sink[predecessor.id] = cost_from_sink[node.id] + cost
                    shortest_path_successor[predecessor.id] = node.id
                    predecessor.shortest_path_successor = node


    def k_shortest_paths(self, k):
        """Generates the k shortest paths from the root to any terminal node.
        
        Generates up to k shortest paths from the root to any terminal node in
        order of increasing length/cost. If less then k paths exist, this
        function returns all paths and stops early.

        Args:
            k (int): The maximum number of paths.

        Yields:
            list: A path, in the form of indices of successors.
            float: The cost of the path.
        """

        # initialization has to happen whenever the graph changed
        if self.k_best_state != len(self.nodes):
            for node in self.nodes:
                path = REAPath(node.optimal_path_predecessor, node.optimal_path_predecessor_idx, 0, node.accumulated_cost)
                node.k_shortest_paths = [path]
                node.k_shortest_paths_candidates = []
            self.k_best_state = len(self.nodes)

        terminal_nodes = list(filter(operator.methodcaller("is_terminal"), self.nodes))
        if not terminal_nodes:
            return
        terminal_nodes.sort(key=operator.attrgetter("accumulated_cost"))

        """
        A graph can have more than one terminal node. It is possible that the
        1st shortest path goes to terminal node 1, the 2nd shortest goes to
        terminal node 2, and the 3rd shortest again goes to terminal node 1.

        The code below ensures that this function generates the shortest paths
        to any terminal node in this graph in order of increasing length,
        irrespective of which terminal node it goes to.

        The idea is to take paths from one terminal node until there is another
        terminal node with a cheaper path. In that case, we continue taking
        paths from that node.
        """

        # contains pairs of [path, generator]
        # this list is sorted by the cost of path because terminal_nodes was sorted
        gens = []
        for node in terminal_nodes:
            gen = node.shortest_paths_iter()
            first_path = next(gen) # The shortest (first) path always exists.
            gens.append([first_path, gen])

        count = 0
        while count < k:
            current_path, gen = gens[0]

            count += 1
            yield current_path

            try:
                next_path = next(gen)
            except StopIteration:
                del gens[0]
                if not gens:
                    return
            else:
                gens[0][0] = next_path

            # if the next path of the current generator is not the cheapest one
            # sort generators
            if len(gens) > 1 and gens[0][0][1] > gens[1][0][1]:
                gens.sort(key=lambda x: x[0][1])

        return


REAPath = collections.namedtuple("REAPath", ["predecessor", "predecessor_idx", "k_prime", "cost"])

class PathDoesNotExist(Exception):
    pass


class GraphNodeBase():

    def __init__(self, factored_operands=None, previous_DS_step=None):

        self.successors = []
        self.edge_labels = []
        self.predecessors = []

        """
        This list holds the indices of self in the successor lists of its
        predecessors. If there are multiple edges to the same predecessor, every
        entry in this list corresponds to a different edge. This is used to
        associate predecessors with edges in case there are multiple edges to
        the same predecessor. More formally:
        For all i, it holds that
        self.predecessors[i].successors[self.predecessors_indices[i]] is self.
        If there is an i and j such that
        self.predecessors[i] is self.predecessors[j], then
        self.predecessors_indices[i] != self.predecessors_indices[j] has to hold.
        """
        self.predecessors_indices = []
        self.original_equations = []
        self.metric = None
        self.accumulated_cost = 0
        self.level = None
        self.optimal_path_predecessor = None
        self.optimal_path_predecessor_idx = None
        self.optimal_path_successor = None
        self.labels = []

        if factored_operands is None:
            self.factored_operands = set()
        else:
            self.factored_operands = factored_operands

        self.applied_DS_steps = set()
        if previous_DS_step:
            self.add_applied_step(previous_DS_step)

        # contains REAPath objects
        self.k_shortest_paths = []
        
        # contains REAPath objects
        self.k_shortest_paths_candidates = []


    def add_factored_operands(self, factored_operands):
        self.factored_operands.update(factored_operands)

    def add_applied_step(self, applied_step):
        """Add an applied derivation step.

        For some derivation steps, the current node stays active. This can cause
        the same derivation step to be applied multiple times to the same node.
        To avoid this, nodes are labelled if those steps are applied.

        Args:
            applied_steps (DS_step): The derivation step that was applied.
        """
        self.applied_DS_steps.add(applied_step)

    def get_payload(self):
        raise NotImplementedError()

    def set_labeled_edge(self, target, label, original_equations):
        idx = len(self.successors)
        self.successors.append(target)
        self.edge_labels.append(label)
        self.original_equations.append(original_equations)

        target.optimal_path_predecessor = self
        target.optimal_path_predecessor_idx = len(target.predecessors)

        # add self to predecessors of target node
        target.predecessors.append(self)
        target.predecessors_indices.append(idx)
        target.accumulated_cost = label.cost + self.accumulated_cost
        # print("set", target.id, self.id, target.accumulated_cost, label.cost, self.accumulated_cost)

    # def remove_edge(self, target):
    #     idx = self.successors.index(target)
    #     self.successors.pop(idx)
    #     self.edge_labels.pop(idx)
    #     self.original_equations.pop(idx)

    def merge(self, other):
        """Merges node other into self.

        After merging, other can (and should) be removed.
        """
        idx_shift = len(self.predecessors)
        for predecessor, p_idx in zip(other.predecessors, other.predecessors_indices):
            predecessor.change_successor(other, self)
            self.predecessors.append(predecessor)
            self.predecessors_indices.append(p_idx)

        for successor in other.successors:
            successor.change_predecessor(other, self)
            if successor.optimal_path_predecessor is other:
                successor.optimal_path_predecessor = self
                # successor.optimal_path_predecessor_idx doesn't have to be changed
        # this can't be done in the loop because change_predecessor uses len(self.successors)
        self.successors.extend(other.successors)

        self.edge_labels.extend(other.edge_labels)
        self.original_equations.extend(other.original_equations)
        self.applied_DS_steps.update(other.applied_DS_steps)
        self.factored_operands.update(other.factored_operands)
        self.update_cost(other, idx_shift)

    def update_cost(self, other, idx_shift):
        """Updates self.accumulated_cost of all successors.

        When merging GraphNodes, the accumulated_cost (which is
        the minimal cost of all paths) of that node most likely changes. If
        there are any successors, those changes have to be propagated
        to all successors. This is what this function does.
        """

        if other.accumulated_cost <= self.accumulated_cost:
            self.accumulated_cost = other.accumulated_cost
            self.optimal_path_predecessor = other.optimal_path_predecessor
            self.optimal_path_predecessor_idx = other.optimal_path_predecessor_idx + idx_shift

            for node in self.topological_sort_successors():
                for idx, (successor, edge_label) in enumerate(zip(node.successors, node.edge_labels)):
                    cost = node.accumulated_cost + edge_label.cost
                    if successor.accumulated_cost > cost:
                        successor.accumulated_cost = cost
                        successor.optimal_path_predecessor = node
                        for new_idx, (predecessor, p_idx) in enumerate(zip(successor.predecessors, successor.predecessors_indices)):
                            if predecessor is node and p_idx == idx:
                                successor.optimal_path_predecessor_idx = new_idx

    def change_successor(self, old_target, new_target):
        idx = self.successors.index(old_target)
        self.successors[idx] = new_target

    def change_predecessor(self, old_target, new_target):
        idx = self.predecessors.index(old_target)
        self.predecessors[idx] = new_target
        self.predecessors_indices[idx] += len(new_target.successors)

    def __str__(self):
        out = "".join(["NODE ", str(self.id), "\n    ", str(self.get_payload()), " ", str(self.metric)])
        predecessor_names = " ".join([predecessor.name for predecessor in self.predecessors])
        out = "".join([out, "\nPREDECESSORS\n    ", predecessor_names, "\nSUCCESSORS"])
        for successor, edge_label in zip(self.successors, self.edge_labels):
            # TODO what is the type of labels?
            out = " ".join([out, "\n    -[", str(edge_label.operation), "]->", str(successor.equations), "(", str(successor.name), ")"])
        return out

    def to_dot(self, style, optimal_edges=set()):
        
        eqns_str = str(self.get_payload())
        # TODO use html module?
        eqns_str = eqns_str.replace('\n', '\\n')
        eqns_str = eqns_str.replace("{", "&#123;")
        eqns_str = eqns_str.replace("}", "&#125;")

        if style is config.GraphStyle.full:
            out = ["""{0} [shape=record, label="{{ {1} |{{ {2} | {3} | {4} | {5:.3g} | {6} | {7} }} }}"];\n""".format(self.name, eqns_str, str(self.id), str(self.level), str(self.metric), self.accumulated_cost, ", ".join([str(op) for op in self.factored_operands]), ", ".join([str(lab) for lab in self.labels]))]
        elif style is config.GraphStyle.simple:
            out = ["""{0} [shape=rect, label="{1}"];""".format(self.name, eqns_str)]
        elif style is config.GraphStyle.minimal:
            out = ["{0} [shape=point];".format(self.name)]

        for successor, label in zip(self.successors, self.edge_labels):
            if style is config.GraphStyle.minimal:
                out.append("{} -> {};".format(self.name, successor.name))
            else:
                if (self.id, successor.id) in optimal_edges:
                    out.append("""{} -> {} [style=bold, label=\"{}\"];\n""".format(self.name, successor.name, str(label)))
                else:
                    out.append("""{} -> {} [label=\"{}\"];\n""".format(self.name, successor.name, str(label)))
        return "".join(out)

    def is_terminal(self):
        raise NotImplementedError()

    def shortest_paths_iter(self):
        """Yields all paths from root to self in order of increasing length.

        Yields:
            list: A path, in the form of indices of successors.
            float: The cost of the path.
        """
        # All paths that have already been computed can be returned right away.
        # This way, paths are not computed more than once, even if this function
        # is called multiple times.
        for path in self.k_shortest_paths:
            yield self.retrieve_path(path)
        # All further paths are computed.
        k = len(self.k_shortest_paths)
        while True:
            try:
                path = self.REA_next_path(k)
            except PathDoesNotExist:
                return
            else:
                k += 1
                yield self.retrieve_path(path)

    def retrieve_path(self, path):
        """Converts REAPath to an actual path.

        Args:
            path (REAPath): The path to convert. It is assumed that it is in
                self.k_shortest_paths.

        Returns:
            list: The converted path, in the form of indices of successors.
            float: The cost of the path.
        """
        _path = []
        current_node = self
        current_path = path
        while current_path.predecessor:
            _path.append(current_node.predecessors_indices[current_path.predecessor_idx])
            current_node = current_path.predecessor
            current_path = current_node.k_shortest_paths[current_path.k_prime]

        return list(reversed(_path)), path.cost

    def REA_next_path(self, k):
        """NextPath(v, k) function of the Recursive Enumeration Algorithm

        This is an implementation of the NextPath(v, k) function of the
        Recursive Enumeration Algorithm for the k shortest paths problem as
        described in the following paper:

        Jiménez, V. M., & Marzal, A. (1999).
        Computing the K Shortest Paths - A New Algorithm and an Experimental
        Comparison. Algorithm Engineering, 1668(Chapter 4), 15–29.
        http://doi.org/10.1007/3-540-48318-7_4

        Quite likely, this implementation has a worse asymptotic complexity than
        the algorithm described in the paper because of the datastructures used
        here and function calls such as index(self) and pop(0). However, so far
        it is sufficiently fast for practical purposes.

        Args:
            self (GraphNodeBase): This is v in the paper.
            k (int): This is k in the paper.

        Returns:
            REAPath: An object that represents the k shortest path.
        """

        # This is not in the paper. It is just for consistency.
        if k == 0:
            return self.k_shortest_paths[0]

        # B.1 initialize candidates
        if k == 1:
            for p_idx, predecessor in enumerate(self.predecessors):
                if p_idx != self.optimal_path_predecessor_idx:
                    idx = self.predecessors_indices[p_idx]
                    cost = predecessor.edge_labels[idx].cost
                    path = REAPath(predecessor, p_idx, 0, predecessor.accumulated_cost + cost)
                    self.k_shortest_paths_candidates.append(path)

        # B.2
        if not (k == 1 and not self.predecessors):
            # B.3
            # get u and k'
            path = self.k_shortest_paths[k-1]
            node_u = path.predecessor
            node_u_idx = path.predecessor_idx
            k_prime = path.k_prime

            # B.4
            if len(node_u.k_shortest_paths) <= k_prime+1:
                try:
                    node_u.REA_next_path(k_prime+1)
                except PathDoesNotExist:
                    pass

            # This can not be the else of the try/except before
            if len(node_u.k_shortest_paths) > k_prime+1:

                # B.5
                # path (k_prime+1, node_u) does exsist

                # get cost
                idx = self.predecessors_indices[node_u_idx]
                cost_uv = node_u.edge_labels[idx].cost
                new_cost = node_u.k_shortest_paths[k_prime+1].cost + cost_uv

                new_path = REAPath(node_u, node_u_idx, k_prime+1, new_cost)
                self.k_shortest_paths_candidates.append(new_path)

        # B.6
        if self.k_shortest_paths_candidates:
            self.k_shortest_paths_candidates.sort(key=operator.attrgetter("cost"))
            path = self.k_shortest_paths_candidates.pop(0)
            if len(self.k_shortest_paths) != k:
                # If this happens, something is completely wrong.
                raise Exception()
            self.k_shortest_paths.append(path)
            return path
        else:
            raise PathDoesNotExist()

    def _topological_sort_visit(self, temp, perm, stack):
        id = self.id
        if id not in temp and id not in perm:
            temp.add(id)
            for successor in self.successors:
                successor._topological_sort_visit(temp, perm, stack)
            perm.add(id)
            temp.remove(id)
            stack.appendleft(self)

    def topological_sort_successors(self):
        stack = collections.deque()
        temp = set()
        perm = set()
        self._topological_sort_visit(temp, perm, stack)
        return list(stack)


class EdgeLabel():
    def __init__(self, *matched_kernels):
        self.matched_kernels = matched_kernels
        self.cost = sum(matched_kernel.cost for matched_kernel in self.matched_kernels)

    def __str__(self):
        lines = []
        for matched_kernel in self.matched_kernels:
            lines.append("{0} {1:.3g}".format(matched_kernel.operation, matched_kernel.cost))
            # lines.append(" ".join([str(matched_kernel.operation), str(int(matched_kernel.cost))]))
        return "\n".join(lines)
