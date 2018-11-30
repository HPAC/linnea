from . import base

from ....algebra.expression import Symbol

class ExpressionGraphBase(base.GraphBase):
    """

    Graph for expressions.
    """

    def __init__(self, root_expr):
        super().__init__()
        self.root = ExpressionGraphNode(root_expr)
        self.active_nodes = [self.root]
        self.nodes = [self.root]


    def create_nodes(self, predecessor, *description):
        new_nodes = []
        # print(description)
        # if description:
        for expression, matched_kernels in description:
            new_node = ExpressionGraphNode(expression, predecessor)
            self.nodes.append(new_node)
            new_nodes.append(new_node)
            predecessor.set_labeled_edge(new_node, base.EdgeLabel(*matched_kernels), None)
        return new_nodes


class ExpressionGraphNode(base.GraphNodeBase):

    _counter = 0

    def __init__(self, expression=None, predecessor=None, factored_operands=None):
        super().__init__(predecessor, factored_operands)
        # IDs for dot output
        self.id = ExpressionGraphNode._counter
        self.name = "".join(["node", str(self.id)])
        ExpressionGraphNode._counter +=1

        self.expression = expression

    def get_payload(self):
        return self.expression

    def is_terminal(self):
        """Returns true if the expression of self is fully computed."""
        return isinstance(self.expression, Symbol)