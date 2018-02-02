
import itertools
import operator
import copy
import math
import os
import collections

from .... import config

output_msg_add = "{2:>2n} {0:.<22}{1:.>5n}"
# output_msg_add = "{2:>2n} {0}:{1:.>27n}"
output_msg_plain = "   {0:.<22}{1:.>5n}"
# output_msg_plain = "   {0}{1:.>27n}"

# TODO rename to GenericGraph or something?
class GraphBase(object):

    def print(self, str):
        if config.verbosity >= 1:
            print(str)

    def print_DS(self, str, val):
        if config.verbosity >= 1:
            print(output_msg_plain.format(str, val))

    def print_DS_numbered(self, str, val, number):
        if config.verbosity >= 1:
            print(output_msg_add.format(str, val, number))

    def derivation(self):
        raise NotImplementedError()

    def nodes_count(self):
        return len(self.nodes)

    def create_nodes(self):
        raise NotImplementedError()

    def remove_nodes(self, nodes):
        """removes nodes from the graph"""
        for node in nodes:
            self.nodes.remove(node)
            try:
                self.active_nodes.remove(node)
            except ValueError:
                pass

    def add_active_nodes(self, nodes):
        self.level_counter += 1
        for node in nodes:
            node.level = self.level_counter

        self.active_nodes.extend(nodes)


    def to_dot(self):
        # TODO use source and sink ranks

        # optimal_nodes = set()
        optimal_edges = set()

        if any(map(operator.methodcaller("is_terminal"), self.nodes)):

            algorithm_paths = list(self.all_algorithms(self.root))
            algorithm_paths.sort(key=operator.itemgetter(1))

            minimal_cost = algorithm_paths[0][1]

            # print(minimal_cost)
            for path, cost in algorithm_paths:
                if cost == minimal_cost:
                    current_node = self.root
                    for path_idx in path:
                        previous_node_id = current_node.id
                        current_node = current_node.successors[path_idx]
                        # optimal_nodes.add(current_node.id)
                        optimal_edges.add((previous_node_id, current_node.id))
                else:
                    break

        # optimal_path = self.optimal_algorithm_path()[0]
        # if optimal_path:
        #     # optimal_nodes.add(self.root.id)
        #     current_node = self.root
        #     for path_idx in optimal_path:
        #         previous_node_id = current_node.id
        #         current_node = current_node.successors[path_idx]
        #         # optimal_nodes.add(current_node.id)
        #         optimal_edges.add((previous_node_id, current_node.id))

        out = "".join([node.to_dot(optimal_edges) for node in self.nodes])
        out = "\n".join(["digraph G {", "ranksep=2.5;", "rankdir=TB;", out, "}"])
        return out


    def to_dot_file(self, name=None):
        if name is "date":
            timestamp = datetime.datetime.now()
            file_name = "".join(["graph_", timestamp.strftime("%Y-%m-%d_%H-%M-%S"), ".gv"])
            output_file = open(file_name, "xt")
        elif name is "counter":
            file_name = "".join(["graph_", str(gn.GraphNode._counter),".gv"])
            output_file = open(file_name, "xt")
        elif name is None:
            file_name = "graph.gv"
            output_file = open(file_name, "wt")
        else:
            file_name = name
            output_file = open(file_name, "wt")
        output_file.write(self.to_dot())
        print("Output was saved in %s" % file_name)
        output_file.close()


    def write_graph(self, output_name, file_name="graph"):
        file_path = os.path.join(config.output_path, output_name, config.language.name, "{}.gv".format(file_name))
        directory_name = os.path.dirname(file_path)
        if not os.path.exists(directory_name):
            os.makedirs(directory_name)
        output_file = open(file_path, "wt")
        output_file.write(self.to_dot())
        output_file.close()
        if config.verbosity >= 2:
            print("Generate graph file {}".format(file_path))

    def __str__(self):
        out = "".join(["Number of nodes: ", str(len(self.nodes)), "\n"])
        out = "".join([out, "Number of active nodes: ", str(len(self.active_nodes)), "\n", "Active nodes: ["])
        nodes = " ".join([node.name for node in self.active_nodes])
        out = "".join([out, nodes, "]"])
        return out

    def all_algorithms(self):
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


    def optimal_algorithm_path(self):
        terminal_nodes = list(filter(operator.methodcaller("is_terminal"), self.nodes))
        if not terminal_nodes:
            return ([], math.inf)

        terminal_nodes.sort(key=operator.attrgetter("accumulated_cost"))
        # print([node.accumulated_cost for node in terminal_nodes])

        current_node = terminal_nodes[0]
        predecessor = current_node.optimal_path_predecessor
        path = []
        while predecessor:
            idx = predecessor.successors.index(current_node)
            path.append(idx)
            current_node = predecessor
            predecessor = current_node.optimal_path_predecessor

        return (list(reversed(path)), terminal_nodes[0].accumulated_cost)


    def optimal_algorithm(self):
        matched_kernels = []
        current_node = self.root
        path, cost = self.optimal_algorithm_path()
        if path:
            for idx in path:

                matched_kernels.extend(current_node.edge_labels[idx].matched_kernels)
                if current_node.successors:
                    current_node = current_node.successors[idx]

            # beginning or end is missing
            return matched_kernels, cost, current_node.get_payload()
        else:
            return [], math.inf, None


    # def _topological_sort_visit(self, node, temp, perm, stack):
    #     if node.id not in temp and node.id not in perm:
    #         node_id = node.id
    #         temp.add(node_id)
    #         for successor in node.successors:
    #             self._topological_sort_visit(successor, temp, perm, stack)
    #         perm.add(node_id)
    #         temp.remove(node_id)
    #         stack.appendleft(node)

    # def topological_sort(self):
    #     stack = collections.deque()
    #     temp = set()
    #     perm = set()
    #     for node in self.nodes:
    #         self._topological_sort_visit(node, temp, perm, stack)

    #     return list(stack)

    # def shortest_path(self):

    #     sorted_nodes = self.topological_sort()


    #     _d = dict()
    #     _d[self.root.id] = 0
    #     # _d = [math.inf]*len(self.nodes)
    #     # _d[0] = 0
    #     _p = dict()
    #     for node in sorted_nodes:
    #         for successor, edge_label in zip(node.successors, node.edge_labels):
    #             cost = edge_label.cost
    #             if _d.get(successor.id, math.inf) > _d[node.id]  + cost:
    #                 _d[successor.id] = _d[node.id] + cost
    #                 _p[successor.id] = node.id

    #     print(_d)
    #     print(_p)
    #     # turns out that _p can easily be obtained by constructing the graph
    #     # in a reasonable way, using gn.GraphNode.optimal_path_predecessor
    #     terminal_nodes = list(filter(is_terminal, self.nodes))

    #     print({node.id: node.optimal_path_predecessor.id for node in self.nodes if node.optimal_path_predecessor})


    # @profile
    def DS_collapse_nodes(self):
        """Merges redundant nodes in the derivation graph.

        Returns the number of removed nodes.
        """

        # This dictionary contains entries of the following structure:
        # 10: [11, 12, 13]
        # This entry means that nodes 10 to 13 are the same.
        duplicates = dict()

        # known_duplicates is needed to avoid {11: [13, 15], 13: [15]}
        known_duplicates = set()

        for p1, p2 in itertools.combinations(enumerate(self.nodes), 2):
            n1, node1 = p1
            n2, node2 = p2
            # The metric is compared first to avoid comparing the actual
            # equations if possible.
            if n2 not in known_duplicates and node1.metric == node2.metric and node1.get_payload() == node2.get_payload():
                duplicates.setdefault(n1, []).append(n2)
                known_duplicates.add(n2)

        #print(duplicates)

        remove = []

        for n in duplicates.keys():
            node = self.nodes[n]
            for i in duplicates[n]:
                redundant_node = self.nodes[i]
                remove.append(redundant_node)
                node.merge(redundant_node)

        self.remove_nodes(remove)

        return len(remove)



class GraphNodeBase(object):

    def __init__(self, metric=None, predecessor=None):

        self.successors = []
        self.edge_labels = []
        self.predecessors = []
        self.metric = metric
        self.accumulated_cost = 0
        self.level = None
        self.optimal_path_predecessor = predecessor

    def get_payload(self):
        raise NotImplementedError()

    def set_labeled_edge(self, target, label):
        self.successors.append(target)
        self.edge_labels.append(label)
        # add self to predecessors of target node
        target.predecessors.append(self)
        target.accumulated_cost = label.cost + self.accumulated_cost
        # print("set", target.id, self.id, target.accumulated_cost, sum(label.cost), self.accumulated_cost)

    def remove_edge(self, target):
        # if target not in self.successors:
        #     # raise exception
        #     pass
        idx = self.successors.index(target)
        self.successors.pop(idx)
        self.edge_labels.pop(idx)

    def merge(self, other):
        """Merges node self with other

        After merging, other can (and should) be removed.
        """
        for predecessor in other.predecessors:
            predecessor.change_successor(other, self)
            self.predecessors.append(predecessor)
        self.successors.extend(other.successors)
        self.edge_labels.extend(other.edge_labels)
        # print("merge", self.id, other.id, self.accumulated_cost, other.accumulated_cost)
        # self.update_cost(min(self.accumulated_cost, other.accumulated_cost))
        self.update_cost(other.accumulated_cost, other.optimal_path_predecessor)
        # self.accumulated_cost = min(self.accumulated_cost, other.accumulated_cost)
        # for successor, edge_label in zip(self.successors, self.edge_labels):
        #     successor.update_cost(sum(edge_label.cost) + self.accumulated_cost)

    def update_cost(self, new_cost, predecessor):
        """Updates self.accumulated_cost of all successors.

        When merging GraphNodes, the accumulated_cost (which is
        the minimal cost of all paths) of that node most likely changes. If
        there are any successors, those changes have to be propagated
        recursively to all successors. This is what this function does.

        """
        if new_cost < self.accumulated_cost:
            self.accumulated_cost = new_cost
            self.optimal_path_predecessor = predecessor
        else:
            return
        for successor, edge_label in zip(self.successors, self.edge_labels):
            successor.update_cost(edge_label.cost + self.accumulated_cost, self)

    def change_successor(self, old_target, new_target):
        idx = self.successors.index(old_target)
        self.successors[idx] = new_target

    def __str__(self):
        out = "".join(["NODE ", str(self.id), "\n    ", str(self.get_payload()), " ", str(self.metric)])
        predecessor_names = " ".join([predecessor.name for predecessor in self.predecessors])
        out = "".join([out, "\nPREDECESSORS\n    ", predecessor_names, "\nSUCCESSORS"])
        for successor, edge_label in zip(self.successors, self.edge_labels):
            # TODO what is the type of labels?
            out = " ".join([out, "\n    -[", str(edge_label.operation), "]->", str(successor.equations), "(", str(successor.name), ")"])
        return out

    def to_dot(self, optimal_edges=set()):
        #out = "".join([self.name, " [shape=box, label=\"", str(self.get_payload()), "\"];\n" ])
        eqns_str = str(self.get_payload())
        # TODO use html module?
        eqns_str = eqns_str.replace('\n', '\\n')
        eqns_str = eqns_str.replace("{", "&#123;")
        eqns_str = eqns_str.replace("}", "&#125;")
        out = """{0} [shape=record, label="{{ {1} |{{ {2} | {3} | {4} | {5:.3g} }} }}"];\n""".format(self.name, eqns_str, str(self.id), str(self.level), str(self.metric), self.accumulated_cost)
        # out = "".join([self.name, " [shape=record, label=\"<f0>", str(self.id), "|<f1>", eqns_str, "|<f2>", str(self.metric), "\"];\n" ])
        for successor, label in zip(self.successors, self.edge_labels):
            if (self.id, successor.id) in optimal_edges:
                # out = "".join([out, self.name, " -> ", successor.name, " [style=bold, color=red, label=\"", str(label), "\"];\n"])
                out = "".join([out, self.name, " -> ", successor.name, " [style=bold, label=\"", str(label), "\"];\n"])
            else:
                out = "".join([out, self.name, " -> ", successor.name, " [label=\"", str(label), "\"];\n"])
        return out

    def is_terminal(self):
        raise NotImplementedError()


class EdgeLabel(object):
    def __init__(self, *matched_kernels):
        self.matched_kernels = matched_kernels
        self.cost = sum(matched_kernel.cost for matched_kernel in self.matched_kernels)

    def __str__(self):
        lines = []
        for matched_kernel in self.matched_kernels:
            lines.append("{0} {1:.3g}".format(matched_kernel.operation, matched_kernel.cost))
            # lines.append(" ".join([str(matched_kernel.operation), str(int(matched_kernel.cost))]))
        return "\n".join(lines)
