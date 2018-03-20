
import matchpy
import itertools

from ..algebra.expression import Times, Plus, Symbol, Operator, \
                                      Transpose, Inverse, InverseTranspose, \
                                      Equal

from ..algebra.equations import Equations

from ..algebra.transformations import invert, invert_transpose, transpose

from ..derivation.graph.utils import is_inverse, is_transpose, is_blocked

from ..algebra.properties import Property as properties

from .. import temporaries

from enum import Enum, unique
from collections import Counter

@unique
class CSEType(Enum):
    none = 0
    inverse = 1
    transpose = 2
    inverse_transpose = 3


def find_max_cliques(G):
    """Returns all maximal cliques in an undirected graph.

    For each node *v*, a *maximal clique for v* is a largest complete
    subgraph containing *v*. The largest maximal clique is sometimes
    called the *maximum clique*.

    This function returns an iterator over cliques, each of which is a
    list of nodes. It is an iterative implementation, so should not
    suffer from recursion depth issues.

    Parameters
    ----------
    G : A dict that represents an undirected graph. They keys are nodes, the
        values are sets of adjacent nodes.

    Returns
    -------
    iterator
        An iterator over maximal cliques, each of which is a list of
        nodes in `G`. The order of cliques is arbitrary.

    Copyright
    --------
    This funciton, including the docstring, is an adapted version of the
    find_cliques function of the NetworkX module.

    Copyright (C) 2004-2018 by
    Aric Hagberg <hagberg@lanl.gov>
    Dan Schult <dschult@colgate.edu>
    Pieter Swart <swart@lanl.gov>
    All rights reserved.
    3-clause BSD license.

    For the full license text, see /other_licenses/networkx.txt

    Notes
    -----
    To obtain a list of all maximal cliques, use
    `list(find_cliques(G))`. However, be aware that in the worst-case,
    the length of this list can be exponential in the number of nodes in
    the graph (for example, when the graph is the complete graph). This
    function avoids storing all cliques in memory by only keeping
    current candidate node lists in memory during its search.

    This implementation is based on the algorithm published by Bron and
    Kerbosch (1973) [1]_, as adapted by Tomita, Tanaka and Takahashi
    (2006) [2]_ and discussed in Cazals and Karande (2008) [3]_. It
    essentially unrolls the recursion used in the references to avoid
    issues of recursion stack depth (for a recursive implementation, see
    :func:`find_cliques_recursive`).

    This algorithm ignores self-loops and parallel edges, since cliques
    are not conventionally defined with such edges.

    References
    ----------
    .. [1] Bron, C. and Kerbosch, J.
       "Algorithm 457: finding all cliques of an undirected graph".
       *Communications of the ACM* 16, 9 (Sep. 1973), 575--577.
       <http://portal.acm.org/citation.cfm?doid=362342.362367>

    .. [2] Etsuji Tomita, Akira Tanaka, Haruhisa Takahashi,
       "The worst-case time complexity for generating all maximal
       cliques and computational experiments",
       *Theoretical Computer Science*, Volume 363, Issue 1,
       Computing and Combinatorics,
       10th Annual International Conference on
       Computing and Combinatorics (COCOON 2004), 25 October 2006, Pages 28--42
       <https://doi.org/10.1016/j.tcs.2006.06.015>

    .. [3] F. Cazals, C. Karande,
       "A note on the problem of reporting maximal cliques",
       *Theoretical Computer Science*,
       Volume 407, Issues 1--3, 6 November 2008, Pages 564--568,
       <https://doi.org/10.1016/j.tcs.2008.05.010>
    """
    if len(G) == 0:
        return

    adj = G
    Q = [None]

    subg = set(G.keys())
    cand = set(G.keys())
    u = max(subg, key=lambda u: len(cand & adj[u]))
    ext_u = cand - adj[u]
    stack = []

    try:
        while True:
            if ext_u:
                q = ext_u.pop()
                cand.remove(q)
                Q[-1] = q
                adj_q = adj[q]
                subg_q = subg & adj_q
                if not subg_q:
                    yield Q[:]
                else:
                    cand_q = cand & adj_q
                    if cand_q:
                        stack.append((subg, cand, ext_u))
                        Q.append(None)
                        subg = subg_q
                        cand = cand_q
                        u = max(subg, key=lambda u: len(cand & adj[u]))
                        ext_u = cand - adj[u]
            else:
                Q.pop()
                subg, cand, ext_u = stack.pop()
    except IndexError:
        pass

class Subexpression(object):
    """docstring for Subexpression"""

    _counter = 0
    def __init__(self, expr, eqn_idx, positions, level, type):
        super(Subexpression, self).__init__()
        self.expr = expr
        self.expr_trans = None
        self.expr_inv = None
        self.expr_invtrans = None
        if type == CSEType.none:
            self.expr_var = [expr]
        elif type == CSEType.inverse:
            self.expr_inv = invert(expr)
            self.expr_var = [self.expr_inv]
        elif type == CSEType.transpose:
            self.expr_trans = transpose(expr)
            self.expr_var = [self.expr_trans]
        elif type == CSEType.inverse_transpose:
            self.expr_inv = invert(expr)
            self.expr_trans = transpose(expr)
            self.expr_invtrans = invert_transpose(expr)
            self.expr_var = [self.expr_inv, self.expr_trans, self.expr_invtrans]

        # print([str(expr) for expr in self.expr_var])
        self.eqn_idx = eqn_idx
        self.positions = positions
        self.type = type
        self.level = level

        self.compatible_count = 0

        self.id = Subexpression._counter
        Subexpression._counter +=1

    def is_compatible(self, other):
        if self.eqn_idx != other.eqn_idx:
            return True
        elif not self.positions.isdisjoint(other.positions):
            return False
        else:
            # check if one position is a prefix of another position
            # if there is more than one position, they all have the same
            # length and only differ in the last position, so it is sufficient
            # to just use one arbitrary position.
            self_pos = next(iter(self.positions))
            other_pos = next(iter(other.positions))

            if len(self_pos) == len(other_pos):
                # We know already that the sets are disjoint. If elements
                # have the same lengths, none is a prefix of the other.
                return True
            else:
                # If one position is prefix of the other, they are not compatible
                # print(self_pos, other_pos)
                for _self_pos, _other_pos in itertools.product(self.positions, other.positions):
                    if all(e1==e2 for e1, e2 in zip(_self_pos, _other_pos)):
                        return False
                return True


    def is_subexpression(self, other):
        """
        This function is transitive.
        """
        if self.eqn_idx != other.eqn_idx:
            return False

        level_diff = self.level - other.level
        if level_diff == 0:
            return self.positions < other.positions
        else:
            # It is sufficient to take an arbitrary positions because they all
            # have the same prefix.
            self_pos = next(iter(self.positions))
            self_prefix = self_pos[0:other.level]
            # print(self_pos, self_prefix, other.positions)
            return self_prefix in other.positions

    def is_transpose(self, other):
        return self.expr == other.expr_trans or self.expr_trans == other.expr
        
    def is_inverse(self, other):
        return self.expr == other.expr_inv or self.expr_inv == other.expr

    def is_inverse_transpose(self, other):
        return self.expr == other.expr_invtrans or self.expr_invtrans == other.expr


def window(seq, n=2):
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    result = tuple(itertools.islice(it, n))
    if len(result) == n:
        yield result    
    for elem in it:
        result = result[1:] + (elem,)
        yield result


def powerset(iterable, min=0, max=None):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    if max is None:
        max = len(s)+1
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(min, max))


# def to_dot_file(relation):
#     out = ""
#     for element in relation:
#         element_str = "{0} [shape=box];\n".format(element.name)
#         out = "".join([out, element_str])

#     for e1, e2 in relation:
#         edge_str = "{0} -> {1};\n".format(self(e2).name, self(e1).name)
#         out = "".join([out, edge_str])        

#     # out = "".join([node.to_dot(optimal_edges) for node in self.nodes])
#     out = "\n".join(["digraph G {", "ranksep=1;", "rankdir=TB;", out, "}"])
    
#     file_name = "partial_ordering.gv"
#     output_file = open(file_name, "wt")
#     output_file.write(out) 
#     print("Output was saved in %s" % file_name)
#     output_file.close()

def contains_inverse(expr):
    if is_inverse(expr):
        return True
    if isinstance(expr, Operator):
        return any(contains_inverse(operand) for operand in expr.operands)

def contains_transpose(expr):
    if is_transpose(expr):
        return True
    if isinstance(expr, Operator):
        return any(contains_transpose(operand) for operand in expr.operands)

def all_subexpressions(expr, position=tuple(), level=0, predeccessor=None):
    """
    Canoical positions: If a subexpression is a single node in the expression
    tree, its position is only one sequence.

    Example: For ABC+D
    ABC+D   {()}                not {(0,), (1,)}
    ABC     {(0,)}              not {(0, 0), (0, 1), (0, 2)}
    D       {(1,)}
    AB      {(0, 0), (0, 1)}

    This way, positions are unique.
    """
    if isinstance(expr, Operator):
        if expr.arity == matchpy.Arity.unary:
            if not (is_inverse(expr) and isinstance(predeccessor, Times)) and not (is_transpose(expr) and isinstance(expr.operand, Symbol)):
                yield expr, {position}, level
        elif expr.commutative:
            if not is_blocked(expr):
                yield expr, {position}, level
            if len(expr.operands) > 2:
                ops = [(op, position + (n,)) for n, op in enumerate(expr.operands)]
                for subset in powerset(ops, 2, len(ops)):
                    positions = set()
                    _ops = []
                    for op, pos in subset:
                        positions.add(pos)
                        _ops.append(op)
                    if not is_blocked(expr):
                        yield Plus(*_ops), positions, level+1
        elif expr.associative:
            if not is_blocked(expr):
                yield expr, {position}, level
            for i in range(2, len(expr.operands)):
                for offset, seq in enumerate(window(expr.operands, i)):
                    positions = set(position + (offset+j,) for j in range(i))
                    if not is_blocked(expr):
                        yield Times(*seq), positions, level+1
        for n, operand in enumerate(expr.operands):
            new_position = position + (n,)
            yield from all_subexpressions(operand, new_position, level+1, expr)



class CSEDetector(object):
    """docstring for CSEDetector"""
    def __init__(self):
        super(CSEDetector, self).__init__()
        self.CSE_detection_dict = dict()
        self.all_subexpressions = dict()
        self.all_CSEs = dict()
        self.subexpr_to_CSE = dict()
        self.all_CSE_lists = []

    def add_subexpression(self, subexpr):

        for expr in subexpr.expr_var:
            try:
                subexprs = self.CSE_detection_dict[expr]
            except KeyError:
                subexprs = []
                self.all_CSE_lists.append(subexprs)
            else:
                break

        subexprs.append(subexpr)
        for expr in subexpr.expr_var:
            self.CSE_detection_dict[expr] = subexprs

        self.all_subexpressions[subexpr.id] = subexpr


    def print_self(self):
        for key, value in self.CSE_detection_dict.items():
            if len(value) > 1:
                print("# ", key)
                for subexpr in value:
                    print(str(subexpr.expr), subexpr.eqn_idx, subexpr.positions)

    def construct_CSEs(self):

        for subexprs in self.all_CSE_lists:
            if len(subexprs) > 1:
                cse = CSE(subexprs)
                self.all_CSEs[cse.id] = cse
                for subexpr in subexprs:
                    self.subexpr_to_CSE[subexpr.id] = cse.id

        print("construct CSEs", [(str(cse.expr), len(cse.subexprs)) for cse in self.all_CSEs.values()])


    def transitive_reduction(self):
        while True:
            remove_edge = []
            for node in self.all_CSEs.values():
                for pre in node.predeccessors:
                    for prepre in pre.predeccessors:
                        if prepre in node.predeccessors:
                            edge = (node, prepre)
                            if edge not in remove_edge:
                                remove_edge.append((node, prepre))
                                # print("tr", str(node.expr), str(prepre.expr))

            if not remove_edge:
                break

            for node1, node2 in remove_edge:
                node1.predeccessors.remove(node2)
                node2.successors.remove(node1)


    def remove_not_maximal_CSEs(self, nodes):
        # TODO do we want to consider compatible_count?

        while True:
            remove = []
            for node in nodes:
                if node.predeccessors and not node.successors:
                    # if not any(len(predeccessor.subexprs) < len(node.subexprs) for predeccessor in node.predeccessors):
                    if not any(predeccessor.max_clique_size < node.max_clique_size for predeccessor in node.predeccessors):
                        # print("not maximal", node.expr)
                        # if not any(predeccessor.compatible_count < node.compatible_count for predeccessor in node.predeccessors):
                        remove.append(node)

            if not remove:
                break

            for node in remove:
                nodes.remove(node)
                for predeccessor in node.predeccessors:
                    predeccessor.successors.remove(node)

    def remove_invalid_CSEs(self, nodes):
        """

        A common subexpression is invalid if it has only two occurrences, and
        those occurrences are not compatible.
        """

        remove = []
        for cse in nodes:
            # if len(cse.subexprs) == 2 and not cse.subexprs[0].is_compatible(cse.subexprs[1]):
            if cse.max_clique_size < 2:
                remove.append(cse)

        for cse in remove:
            nodes.remove(cse)
            for predeccessor in cse.predeccessors:
                predeccessor.successors.remove(cse)
            for successor in cse.successors:
                successor.predeccessors.remove(cse)


    def maximal_CSEs(self):
        current_nodes = [node for node in self.all_CSEs.values()]

        # print("all nodes", [str(node.expr) for node in current_nodes])
        # remove not maximal CSEs

        self.remove_invalid_CSEs(current_nodes)

        self.remove_not_maximal_CSEs(current_nodes)

        # print("all nodes", [str(node.expr) for node in current_nodes])

        CSEs = []
        while current_nodes:
            # print("lattice", [str(node.expr) for node in lattice])
            new_CSEs = []
            for node in current_nodes:
                # print("node", str(node.expr))
                # print("successors", [str(node.expr) for node in node.successors])
                if not node.successors:
                    new_CSEs.append(node)

            # print("new_CSEs", [str(node.expr) for node in new_CSEs])
            for node in new_CSEs:
                current_nodes.remove(node)
                for predeccessor in node.predeccessors:
                    predeccessor.successors.remove(node)

            self.remove_not_maximal_CSEs(current_nodes)

            CSEs.extend(new_CSEs)

        return CSEs


    def replaceable_CSEs(self, CSEs):

        if len(CSEs) == 1:
            # No constraints need to be checked here. All subexpressions belong
            # to the same CSE, and there are no cliques with less than two
            # nodes.
            for clique in CSEs[0].all_cliques:
                yield [self.all_subexpressions[id] for id in clique]
        else:
            subexprs = list(itertools.chain.from_iterable([node.subexprs for node in CSEs]))

            adj = dict()
            for subexpr1, subexpr2 in itertools.combinations(subexprs, 2):
                if subexpr1.is_compatible(subexpr2):
                    # print(subexpr1.expr, subexpr2.expr)
                    adj.setdefault(subexpr1.id, set()).add(subexpr2.id)
                    adj.setdefault(subexpr2.id, set()).add(subexpr1.id)

            cliques = []
            for clique in find_max_cliques(adj):

                _counter = Counter(self.subexpr_to_CSE[id] for id in clique)
                # print(clique)
                # print(_counter)

                # if for a CSE there is only one subexpression, this subexpression is removed
                remove = set()
                for id in clique:
                    if _counter[self.subexpr_to_CSE[id]] == 1:
                        remove.add(id)
                clique_no_singles = set(clique) - remove
                cliques.append(clique_no_singles)

            # removing cliques that are subsets of other cliques
            # this is possible because of the previous step
            remove = []
            for clique1, clique2 in itertools.combinations(cliques, 2):
                if clique1 <= clique2:
                    remove.append(clique1)
                elif clique2 < clique1:
                    remove.append(clique2)

            for clique in remove:
                try:
                    cliques.remove(clique)
                except ValueError:
                    pass
                
            # print(cliques)

            for clique in cliques:
                yield [self.all_subexpressions[id] for id in clique]


    def construct_lattice(self):
        for cse1, cse2 in itertools.product(self.all_CSEs.values(), repeat=2):
            if cse1.is_subexpression(cse2):
                cse1.predeccessors.append(cse2)
                cse2.successors.append(cse1)
                # print(str(cse1.expr), str(cse2.expr))

    def count_compatible(self):
        
        # print("counting other")
        for cse1, cse2 in itertools.combinations(self.all_CSEs.values(), 2):
            # print("CSE", str(cse1.expr), str(cse2.expr), len(cse2.subexprs))
            for subexpr1, subexpr2 in itertools.product(cse1.subexprs, cse2.subexprs):
                if subexpr1.is_compatible(subexpr2):
                    # print(str(subexpr1.expr), subexpr1.eqn_idx, subexpr1.positions, str(subexpr2.expr), subexpr2.eqn_idx, subexpr2.positions)
                    subexpr1.compatible_count += 1
                    subexpr2.compatible_count += 1

        # print("counting self")
        for cse in self.all_CSEs.values():
            for subexpr1, subexpr2 in itertools.combinations(cse.subexprs, 2):
                if subexpr1.is_compatible(subexpr2):
                    # print(str(subexpr1.expr), subexpr1.eqn_idx, str(subexpr2.expr), subexpr2.eqn_idx)
                    subexpr1.compatible_count += 1
                    subexpr2.compatible_count += 1


        for CSE in self.all_CSEs.values():
            CSE.compatible_count = sum(subexpr.compatible_count for subexpr in CSE.subexprs)

        print("compatible count", [(str(cse.expr), cse.compatible_count) for cse in self.all_CSEs.values()])


    def CSEs(self):
        self.construct_CSEs()

        self.construct_lattice()

        # self.count_compatible()

        self.transitive_reduction()
        maximal_CSEs = self.maximal_CSEs()

        # adding pairs of CSEs that are compatible
        # there is no reason not to go for larger groups, but they are probably very unlikely
        grouped_CSEs = []
        for cse1, cse2 in itertools.combinations(self.all_CSEs.values(), 2):
            if cse1.is_compatible(cse2):
                grouped_CSEs.append((cse1, cse2))

        for group in grouped_CSEs:
            for cse in group:
                try:
                    maximal_CSEs.remove(cse)
                except ValueError:
                    pass

        for cse in maximal_CSEs:
            grouped_CSEs.append((cse,))      

        for group in grouped_CSEs:
            print("group", [str(cse.expr) for cse in group])
            yield from self.replaceable_CSEs(group)



class CSE(object):
    """docstring for CSE"""

    _counter = 0
    def __init__(self, subexprs):
        super(CSE, self).__init__()
        self.subexprs = subexprs
        self.subexprs_ids = {subexpr.id for subexpr in self.subexprs}
        self.expr = subexprs[0].expr # for DEBUG only

        self.id = CSE._counter
        CSE._counter +=1

        self.compatible_count = 0

        self.predeccessors = []
        self.successors = []

        adj = dict()
        for subexpr1, subexpr2 in itertools.combinations(self.subexprs, 2):
            if subexpr1.is_compatible(subexpr2):
                # print(subexpr1.expr, subexpr2.expr)
                adj.setdefault(subexpr1.id, set()).add(subexpr2.id)
                adj.setdefault(subexpr2.id, set()).add(subexpr1.id)

        self.all_cliques = list(find_max_cliques(adj))
        self.max_clique_size = max((len(clique) for clique in self.all_cliques), default=0)

    def is_subexpression(self, other):
        for self_subexpr, other_subexpr in itertools.product(self.subexprs, other.subexprs):
            if self_subexpr.is_subexpression(other_subexpr):
                return True
        return False

    def is_compatible(self, other):
        for self_subexpr, other_subexpr in itertools.product(self.subexprs, other.subexprs):
            if not self_subexpr.is_compatible(other_subexpr):
                return False
        return True


def indentify_subexpression_types(subexprs):
    # subexprs_dict = dict()
    # for subexpr in subexprs:
    #     subexprs_dict.setdefault(subexpr.expr, []).append(subexpr)

    # # print(subexprs_dict)

    # if len(subexprs_dict.keys()) == 1:
    #     return [CSEType.none for _ in subexprs]
    # else:

    ref_subexpr = subexprs[0]
    return_types = [CSEType.none]
    for subexpr in subexprs[1:]:
        if subexpr.expr == ref_subexpr.expr:
            return_types.append(CSEType.none)
        elif subexpr.is_transpose(ref_subexpr):
            return_types.append(CSEType.transpose)
        elif subexpr.is_inverse(ref_subexpr):
            return_types.append(CSEType.inverse)
        elif subexpr.is_inverse_transpose(ref_subexpr):
            return_types.append(CSEType.invert_transpose)
    
    # print(return_types)
    return return_types


def find_CSEs(equations):


    CSE_detector = CSEDetector()
    for eqn_idx, equation in enumerate(equations):
        for expr, positions, level in all_subexpressions(equation.rhs):
            # print(expr, positions)

            # for expression of the form Inverse(expr), we don't want to add the
            # inverted variant because this will produce "fake" CSEs. The reason
            # is that expr will also be added.
            # same for transpose
            # TODO do we want to invert if expr is not square?
            inv = contains_inverse(expr) and not is_inverse(expr)
            trans = contains_transpose(expr) and not is_transpose(expr)
            if expr.has_property(properties.SYMMETRIC):
                trans = False
            if inv and trans:
                subexpr = Subexpression(expr, eqn_idx, positions, level, CSEType.inverse_transpose)
                # print("inv trans", expr)
                CSE_detector.add_subexpression(subexpr)
            elif inv:
                subexpr = Subexpression(expr, eqn_idx, positions, level, CSEType.inverse)
                # print("inv", expr, subexpr.expr_var)
                CSE_detector.add_subexpression(subexpr)
            elif trans:
                subexpr = Subexpression(expr, eqn_idx, positions, level, CSEType.transpose)
                # print("trans", expr, subexpr.expr_var)
                CSE_detector.add_subexpression(subexpr)
            else:
                subexpr = Subexpression(expr, eqn_idx, positions, level, CSEType.none)
                CSE_detector.add_subexpression(subexpr)


    # CSE_detector.print_self()
    for CSE in CSE_detector.CSEs():
        # pass
        print("CSEs", [(str(subexpr.expr), subexpr.eqn_idx) for subexpr in CSE])

        CSE_as_dict = dict()
        for subexpr in CSE:
            CSE_as_dict.setdefault(CSE_detector.subexpr_to_CSE[subexpr.id], []).append(subexpr)
  
        insert_equations = []
        replacements_per_equation = dict()
        for CSE_id, subexprs in CSE_as_dict.items():

            min_eqn_idx = min(subexpr.eqn_idx for subexpr in subexprs)

            CSE_expr = subexprs[0].expr # this works because indentify_subexpression_types uses subexprs[0] as reference
            tmp = temporaries.create_tmp(CSE_expr, True)
            eqn = Equal(tmp, CSE_expr)

            insert_equations.append((min_eqn_idx, eqn))

            for subexpr, subexpr_type in zip(subexprs, indentify_subexpression_types(subexprs)):
                positions = list(subexpr.positions)
                for position in positions[1:]:
                    replacements_per_equation.setdefault(subexpr.eqn_idx, []).append(((1,) + position, []))

                tmp_expr = tmp
                if subexpr_type == CSEType.transpose:
                    tmp_expr = Transpose(tmp)
                elif subexpr_type == CSEType.inverse:
                    tmp_expr = Inverse(tmp)
                elif subexpr_type == CSEType.inverse_transpose:
                    tmp_expr = InverseTranspose(tmp)

                replacements_per_equation.setdefault(subexpr.eqn_idx, []).append(((1,) + positions[0], tmp_expr))


        # print(replacements_per_equation)
        equations_list = list(equations.equations)
        for eqn_idx, replacements in replacements_per_equation.items():
            equations_list[eqn_idx] = matchpy.replace_many(equations_list[eqn_idx], replacements)
            

        # Inserting new equation for extracted CSE.
        # It is inserted right before the first occurrence of the CSE.
        insert_equations.sort(reverse=True)
        for min_eqn_idx, eqn in insert_equations:
            equations_list.insert(min_eqn_idx, eqn)
        new_equations = Equations(*equations_list)
        new_equations = new_equations.to_normalform()

        print(new_equations)

        yield new_equations



if __name__ == "__main__":

    import linnea.config

    linnea.config.init()

    from linnea.derivation.graph.constructive import DerivationGraph
    # from linnea.derivation.graph.exhaustive import DerivationGraph
    # from linnea.derivation.graph.matrix_chain_derivation import DerivationGraph

    import linnea.examples.examples
    import linnea.examples.lamp_paper.examples

    # example = linnea.examples.lamp_paper.examples.Common_Subexpr_7_2_3()
    # example = linnea.examples.lamp_paper.examples.Overlap_Common_Subexpr_7_2_4()
    # example = linnea.examples.lamp_paper.examples.Matrix_Chain_7_2_6()
    # example = linnea.examples.lamp_paper.examples.Generalized_LeastSquares_7_1_3()
    # example = linnea.examples.lamp_paper.examples.Rewrite_Distributivity_7_2_5_1()
    # example = linnea.examples.lamp_paper.examples.LMMSE_7_1_2() # TODO exhaustive: eigen trick is not applied?
    # example = linnea.examples.lamp_paper.examples.LeastSquares_7_1_1()
    # example = linnea.examples.lamp_paper.examples.Simplification_7_2_12()
    # example = linnea.examples.lamp_paper.examples.Simplification_7_2_11()
    # example = linnea.examples.lamp_paper.examples.Tikhonov_7_1_14() # TODO exhaustive: eigen trick is not applied because of dead ends?
    # example = linnea.examples.lamp_paper.examples.EnsembleKalmarFilter_7_1_9_1() # large number of solutions
    # example = linnea.examples.lamp_paper.examples.EnsembleKalmarFilter_7_1_9_2()
    # example = linnea.examples.lamp_paper.examples.Common_Subexpr_7_2_2()
    # example = linnea.examples.lamp_paper.examples.CDMA_7_1_15() # set_equivalent_upwards makes a big difference for number of nodes here
    # example = linnea.examples.lamp_paper.examples.Local_Assimilation_Kalmar_7_1_7() # TODO problem with storage format (constructive)
    # example = linnea.examples.lamp_paper.examples.ImageRestoration_7_1_13_1() # TODO POS problem? CSE problem
    # example = linnea.examples.lamp_paper.examples.ImageRestoration_7_1_13_2() 
    # example = linnea.examples.lamp_paper.examples.ImageRestoration_7_1_13_2(single=True) # all_algorithms is expensive here with exhaustive and dead_ends=False

    # CSE example

    # example = linnea.examples.examples.Example002()
    # example = linnea.examples.examples.Example004()
    # example = linnea.examples.examples.Example009()
    # example = linnea.examples.examples.Example014()
    # example = linnea.examples.examples.Example025()
    # example = linnea.examples.examples.Example028()
    # example = linnea.examples.examples.Example032()
    # example = linnea.examples.examples.Example033()
    # example = linnea.examples.examples.Example040()
    # example = linnea.examples.examples.Example041()
    # example = linnea.examples.examples.Example043()
    # example = linnea.examples.examples.Example044()
    # example = linnea.examples.examples.Example045()
    # example = linnea.examples.examples.Example046()
    # example = linnea.examples.examples.Example050()
    # example = linnea.examples.examples.Example051()
    # example = linnea.examples.examples.Example052()
    # example = linnea.examples.examples.Example053()
    # example = linnea.examples.examples.Example055() # TODO do we want C^-1 B^-1 ?
    # example = linnea.examples.examples.Example056() # Plus
    # example = linnea.examples.examples.Example057() # count compatible CSEs
    # example = linnea.examples.examples.Example058()
    # example = linnea.examples.examples.Example065()
    example = linnea.examples.examples.Example066() #!! maximal CSE
    # example = linnea.examples.examples.Example067()
    example = linnea.examples.examples.Example069() # count compatible works (because of counting with itself)
    # example = linnea.examples.examples.Example078()
    # example = linnea.examples.examples.Example079()
    # example = linnea.examples.examples.Example080() # type
    # example = linnea.examples.examples.Example081() # type TODO do we want B^-T C^-T
    # example = linnea.examples.examples.Example082() # type
    # example = linnea.examples.examples.Example083() # type
    # example = linnea.examples.examples.Example085() # count compatible works (because of counting with itself)
    # example = linnea.examples.examples.Example086() # count compatible works (because of counting with itself)
    # example = linnea.examples.examples.Example093() # overlapping CSEs # TODO do we want A B^-1 here?
    # example = linnea.examples.examples.Example104() # explicit inversion
    # example = linnea.examples.examples.Example116()
    # example = linnea.examples.examples.Example124() # type
    # example = linnea.examples.examples.Example130()
    # example = linnea.examples.examples.Example131()
    # example = linnea.examples.examples.Example132()
    # example = linnea.examples.examples.Example134() # b = ((A (X^T X)^-1 X^T y) + (B (X^T X)^-1 X^T z))
    # example = linnea.examples.examples.Example135() #!! on paper
    # example = linnea.examples.examples.Example136() #!! even harder: removing wrong cliques here is difficult -> subgraphs
    # example = linnea.examples.examples.Example137() #!! maximal CSE
    # example = linnea.examples.examples.Example138() #!! maximal CSE # count compatible works
    # example = linnea.examples.examples.Example139() # count compatible works
    # example = linnea.examples.examples.Example140() # count compatible works
    # example = linnea.examples.examples.Example141() # lots of solutions, count compatible, maximal CSE, connected components
    # example = linnea.examples.examples.Example142() # overlapping with itself

    eqns = example.eqns.to_normalform()

    print(eqns)

    for res in find_CSEs(eqns):
        pass

    # for equation in example.eqns:
    #     # print(equation.rhs)
    #     for subexpr, pos in all_subexpressions(equation.rhs):
    #         print(str(subexpr), pos)


    # graph = networkx.Graph()
    # graph.add_edges_from({(0,1), (1,2), (2,3),(3,1)})
    # # print(graph.nodes())
    # for clique in networkx.find_cliques(graph):
    #     print(clique)
