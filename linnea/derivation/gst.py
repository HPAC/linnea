from collections import deque

class GST_node(object):

    _counter = 0

    def __init__(self):
        
        self.id = GST_node._counter
        self.name = "".join(["node", str(self.id)])
        GST_node._counter += 1

        self.successors = []
        self.edge_labels = []
        self.suffix_link = None
        self.reverse_suffix_links = []
        # Used to identify all nodes on the path of the full input sequence
        self.full_path = False

        self.positions = dict()

    def set_successor(self, successor, edge_label):
        self.successors.append(successor)
        self.edge_labels.append(edge_label)

    def set_suffix_link(self, target):
        """Sets suffix link of self to target (plus reversed link)"""
        self.suffix_link = target
        target.reverse_suffix_links.append(self)

    def merge_positions(self, other):
        """Merges self.positions and other.positions"""
        for key in other.positions.keys():
            self.positions.setdefault(key, []).extend(other.positions[key])

    def get_successor(self, symbol):
        """Returns the successor of self following the edge labelled with symbol."""
        try:
            idx = self.edge_labels.index(symbol)
        except ValueError:
            return None
        else:
            return self.successors[idx]

    def to_dot(self):
        if self.full_path:
            style = "style=filled, color=lightgray,"
        else:
            style = ""
        out = "".join([self.name, " [shape=box, ", style, " label=\"", str(self.id), "\n", str(self.positions), "\"];\n" ])
        # out = "".join([self.name, " [shape=box, ", style, " label=\"", str(self.indices), "\n", str(self.id), "\"];\n" ])
        for successor, label in zip(self.successors, self.edge_labels):
            out = "".join([out, self.name, " -> ", successor.name, " [label=\"", str(label), "\"];\n"])
        
        # if self.suffix_link:
        #     out = "".join([out, self.name, " -> ", self.suffix_link.name, " [style=dotted];\n"])
        if self.reverse_suffix_links:
            for link in self.reverse_suffix_links:
                out = "".join([out, self.name, " -> ", link.name, " [style=dashed];\n"])
        return out

class GST(object):
    def __init__(self, sequences=None):
        self.root = GST_node()
        self.root.full_path = True
        self.nodes = [self.root]

        # Counts the input sequences.
        self.counter = 0

        # Stores the length of each input sequence.
        self.input_lengths = []

        # This makes sure that positions are only updated once.
        # TODO remove?
        # There is no reason to compue CSEs more than once.
        self.update_positions = True

        if sequences:
            for sequence in sequences:
                self.add_sequence(sequence)

    def add_sequences(self, sequences):
        for sequence in sequences:
            self.add_sequence(sequence)

    

    # def longest_replaceable_CSE(self):
    #     """

    #     A CSE is replaceable if there are at least two occurences that don't
    #     overlap with itself. Example:
    #     AAAA
    #     The CSE "AAA" is not replaceable, because the occurences at positions 0
    #     and 1 overlap.
    #     """
    #     CSEs = self.find_all_CSEs()
    #     CSEs.sort(key=lambda x: x[1], reverse=True)

        
    #     for positions, length in CSEs:
    #         remaining_positions = dict()    
    #         for seq_idx, starting_pos in positions.items():
    #             if len(starting_pos) > 1:
    #                 starting_pos.sort()
    #                 for n, pos in enumerate(starting_pos[:-1]):
    #                     if pos + length <= starting_pos[n+1]:
    #                         remaining_positions.setdefault(seq_idx, []).append(pos)
    #             else:
    #                 remaining_positions[seq_idx] = starting_pos
    #         print(remaining_positions)

        # for positions, length in CSEs:
        #     for seq_idx in positions.keys():
        #         starting_pos = positions[seq_idx]
        #         if len(starting_pos) > 1:
        #             starting_pos.sort()
        #             for n, pos in enumerate(starting_pos[:-1]):
        #                 if pos + length < starting_pos[n+1]:

        #             print(starting_pos)


    def find_all_CSEs(self, min_length=2):
        """ Finds all common subexpression with length >= min_length.

        Returns a list of tuples. Each tuple contains
        - a dictionary
        - the length of this CSE.

        The dictionary is used as follows:
        {0: [1, 4], 1: [3]}
        means that the CSE appears in the first input sequence (0) at position 1
        and 4 (starting at 0) and in the second input sequence (1) at position 3.
        """

        CSEs = []

        for node, length in self._find_all_CSEs(self.root, 0, min_length):
            CSEs.append((node.positions, length))

        self.update_positions = False

        return CSEs

    def _find_all_CSEs(self, node, depth, min_length):
        """Finds all nodes in the tree that represents CSEs.

        Yields pairs with
        - a node that represents a CSE
        - the depth at which that node was found
        """
        # Traversing the entire tree.
        for successor in node.successors:
            yield from self._find_all_CSEs(successor, depth+1, min_length)
            if self.update_positions:
                # Updating positions along the way.
                node.merge_positions(successor)

        # A node
        # - with two or more incoming suffix links OR
        # - that lies on the path of the full input sequence and has one incoming
        #   suffix link
        # - that has more than one descendant where a suffix ends (> 1 entries in
        #   node.positions) # TODO does that make 1 and 2 redundant?
        #   Remark: With this, we also find CSEs that are prefixes of longer CSEs.
        # - COMMENTED OUT that is a leaf where more than one suffix ends (> 1 entries in
        #   node.positions)
        # represents a common subexpression.

        reversed_links = len(node.reverse_suffix_links)

        number_of_position = sum(len(positions) for positions in node.positions.values())

        # I'm not sure if this is correct, but it seems to work just fine.
        if depth >= min_length:
            if number_of_position > 1:
                yield (node, depth)

        # if (reversed_links >= 2 or
        #    (reversed_links == 1 and node.full_path)
        #    or (len(node.positions) > 1)
        #    ) and depth >= min_length:
        #     yield (node, depth)

        # In theory, this selects only maximal CSEs. Maximal means that if a CSE
        # is extended, then one or more occurrences are lost.
        # See "Structural Analysis of Gapped Motifs of a String" by Esko
        # Ukkonen, section "Representation by Maximal Non–gapped Motifs"
        # The problem is that in case of expr = A^T B^T B A, the CSE B A is not
        # maximal if both expr and expr^T are in the tree.
        # Looks like the only way to avoid non-maximal CSEs is to sort them out
        # later.
        # TODO Using len(node.positions) is not correct. If there are mutliple
        # occurrences in the same sequence, this doesn't count all occurrences
        # because they are stored as a dictionary where the values are lists and
        # each entry in one list is one occurrence.
        # if depth >= min_length and len(node.positions) > 1:
        #     if node.full_path:
        #         if len(node.successors) == 1:
        #             if len(node.positions) > len(node.successors[0].positions):
        #                 yield (node, depth)
        #         else: # this means len(node.successors) is 0 or larger 1. Is 0 okay?
        #             yield (node, depth)
        #     else:
        #         if reversed_links > 1:
        #             yield (node, depth)
    
    def add_sequence(self, sequence):
        """This function implements Ukkonen's algorithm."""

        self.update_positions = True

        # Don't modify input
        sequence = sequence.copy()
        n = len(sequence)
        self.input_lengths.append(n)
        deepest_node = self.root
        while sequence:
            symbol = sequence.pop(0)

            # If there already is a node for node for the current symbol, go to
            # that node and skip the rest. This happens only if there is more 
            # than one sequence
            successor = deepest_node.get_successor(symbol)
            if successor:
                successor.full_path = True
                deepest_node = successor
                continue

            new_node1 = GST_node()
            new_node1.full_path = True
            self.nodes.append(new_node1)
            deepest_node.set_successor(new_node1, symbol)

            current_node = deepest_node
            deepest_node = new_node1

            while True:
                # This is the case only if current_node is root.
                if not current_node.suffix_link:
                    new_node1.set_suffix_link(current_node)
                    break
                # Following the suffix links.
                current_node = current_node.suffix_link
                
                if symbol not in current_node.edge_labels:
                    # A new node has to be created
                    new_node2 = GST_node()
                    self.nodes.append(new_node2)
                    current_node.set_successor(new_node2, symbol)
                    new_node1.set_suffix_link(new_node2)
                else:
                    # If a node exist already, we set the suffix link and stop.
                    # All other suffix links allong the suffix link chain are 
                    # already set.
                    successor = current_node.get_successor(symbol)
                    new_node1.set_suffix_link(successor)
                    break
                new_node1 = new_node2

        # Here, we traverse the suffix link chain and mark all nodes along it
        # as terminal nodes for the current input. Clearly, this happens only
        # once per input sequence, when it is fully processed (yes, this is
        # correct).
        k = 0
        current_node = deepest_node
        while current_node != self.root:
            current_node.positions.setdefault(self.counter, []).append(k)
            k += 1
            current_node = current_node.suffix_link
        self.counter += 1

    def to_dot(self):
        out = "".join([node.to_dot() for node in self.nodes])
        out = "\n".join(["digraph G {", "ranksep=1;", "rankdir=TB;", out, "}"])
        return out

    def to_dot_file(self):
        file_name = "gst.gv"
        output_file = open(file_name, "wt")
        output_file.write(self.to_dot())
        print("Output was saved in %s" % file_name)