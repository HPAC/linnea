from . import base

from ....algebra.expression import Symbol

class ExpressionGraphBase(base.GraphBase):
    """

    Graph for expressions.
    """

    def __init__(self, root_expr):
        super(ExpressionGraphBase, self).__init__()
        self.root = ExpressionGraphNode(root_expr)
        self.active_nodes = [self.root]
        self.nodes = [self.root]
        self.level_counter = -1


    def create_nodes(self, predecessor, *description):
        new_nodes = []
        # print(description)
        # if description:
        for expression, metric, edge_label in description:
            new_node = ExpressionGraphNode(expression, metric, predecessor)
            self.nodes.append(new_node)
            new_nodes.append(new_node)
            predecessor.set_labeled_edge(new_node, edge_label)
        return new_nodes


class ExpressionGraphNode(base.GraphNodeBase):

    _counter = 0

    def __init__(self, expression=None, metric=None, predecessor=None):
        super(ExpressionGraphNode, self).__init__(metric, predecessor)
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