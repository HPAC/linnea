from ..algebra.expression import Times, Plus, Symbol, Operator, \
                                 Transpose, Inverse, InverseTranspose, \
                                 Equal
from ..algebra.equations import Equations
from ..algebra.transformations import invert, invert_transpose, transpose
from ..algebra.properties import Property
from ..utils import window, powerset, is_inverse, is_transpose, \
                    contains_inverse, contains_transpose

from .. import temporaries

from .utils import is_blocked

from enum import Enum, unique
from collections import Counter

import matchpy
import itertools

@unique
class CSEType(Enum):
    none = 0
    inverse = 1
    transpose = 2
    inverse_transpose = 3

class Subexpression():
    """Subexpression for CSE detection.

    This class represents a subexpression to be used for CSE detection.
    """

    _counter = 0
    def __init__(self, expr, eqn_idx, positions, level, type):
        super().__init__()
        self.expr = expr
        self.expr_trans = None
        self.expr_inv = None
        self.expr_invtrans = None
        if type == CSEType.none:
            self.expr_var = [self.expr]
        elif type == CSEType.inverse:
            self.expr_inv = invert(expr)
            self.expr_var = [self.expr, self.expr_inv]
        elif type == CSEType.transpose:
            self.expr_trans = transpose(expr)
            self.expr_var = [self.expr, self.expr_trans]
        elif type == CSEType.inverse_transpose:
            self.expr_inv = invert(expr)
            self.expr_trans = transpose(expr)
            self.expr_invtrans = invert_transpose(expr)
            self.expr_var = [self.expr, self.expr_inv, self.expr_trans, self.expr_invtrans]

        # print([str(expr) for expr in self.expr_var])
        self.eqn_idx = eqn_idx
        self.positions = positions
        self.type = type
        self.level = level

        self.id = Subexpression._counter
        Subexpression._counter +=1

    def __lt__(self, other):
        return self.id < other.id

    def is_compatible(self, other):
        """Tests if self is compatible with other.

        Two subexpressions are compatible if they do not overlap, or if none is
        a subexpression of the other. This means that both subexpressions can
        be replaced at the same time.

        Args:
            other (Subexpression): The other subexpression.

        Returns:
            True if self is compatible with other, False otherwise.       
        """
        if self.eqn_idx != other.eqn_idx:
            return True
        elif not self.positions.isdisjoint(other.positions):
            return False
        else:
            # check if one position is a prefix of another position
            # if there is more than one position, they all have the same
            # length and only differ in the last position, so it is sufficient
            # to just use one arbitrary position.

            if self.level == other.level:
                # We know already that the sets are disjoint. If elements
                # have the same lengths, none is a prefix of the other.
                return True
            else:
                # If one position is prefix of the other, they are not compatible
                # print(self_pos, other_pos)
                # TODO it should be possible to do this without product.
                # Check which one is longer, then check prefix.
                for self_pos, other_pos in itertools.product(self.positions, other.positions):
                    """This is a simple way to check the prefix relatioship.
                    zip stops as soon as the first iterable is exhausted, so we
                    are only looking at the elements that both may have in
                    common. The case where one path (i.e. iterable) is empty is
                    handled correctly because all returns true for an empty
                    iterable.
                    """
                    if all(e1==e2 for e1, e2 in zip(self_pos, other_pos)):
                        return False
                return True


    def is_subexpression(self, other):
        """Tests if self is a subexpression of other.

        The test is based on positions, not expressions, so variants of
        expressions do not play a role. However, his also means that this
        function returns False if the self.expr is an actual subexpression of
        other.expr in case they come from different positions.

        This function defines a transitive relation.

        Args:
            other (Subexpression): The other subexpression.

        Returns:
            True if self is a subexpression of other, False otherwise.
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
        """Tests if self is the transpose of other.

        Args:
            other (Subexpression): The other subexpression.

        Returns:
            True if self is the transpose of other, False otherwise.      
        """
        return self.expr == other.expr_trans or self.expr_trans == other.expr
        
    def is_inverse(self, other):
        """Tests if self is the inverse of other.

        Args:
            other (Subexpression): The other subexpression.

        Returns:
            True if self is the inverse of other, False otherwise.      
        """
        return self.expr == other.expr_inv or self.expr_inv == other.expr

    def is_inverse_transpose(self, other):
        """Tests if self is the inverse-transpose of other.

        Args:
            other (Subexpression): The other subexpression.

        Returns:
            True if self is the inverse-transpose of other, False otherwise.      
        """
        return self.expr == other.expr_invtrans or self.expr_invtrans == other.expr

class CSE():
    """Common Subexpression object

    This is mostly just a collection of Subexpression objects.
    """

    _counter = 0
    def __init__(self, subexprs):
        super().__init__()
        self.subexprs = subexprs
        self.subexprs_ids = {subexpr.id for subexpr in self.subexprs}
        self.expr = subexprs[0].expr # for DEBUG only

        self.id = CSE._counter
        CSE._counter +=1

        adj = construct_graph(self.subexprs)
        self.all_cliques = list(find_max_cliques(adj))
        """int: This number gives the maximal number of subexpressions that can
        be replaced at the same time. It is given by the size the largest
        maximal clique in a graph where subexpressions are nodes, and edges
        represents which subexpressions are compatible.
        I suspect that for products, those graphs are (proper) intervals graphs.
        For interval graphs, cliques can be computed in polynomial time. Do
        graphs for sums have special properties too?
        """
        self.max_clique_size = max((len(clique) for clique in self.all_cliques), default=0)

    def is_subexpression(self, other):
        """Tests if self is a subexpression of other.

        The test is based on the positions (occurrences) of the subexpressions
        of self and other. If for every occurrence in self, there is an
        occurrence in other where the occurrence in self is subexpression of the
        occurrence in other, this function returns True.

        This function defines a transitive relation.

        Args:
            other (CSE): The other CSE

        Returns:
            True if self is a subexpression of other, False otherwise.
        """
        return all(any(self_subexpr.is_subexpression(other_subexpr) for other_subexpr in other.subexprs) for self_subexpr in self.subexprs)

    def is_compatible(self, other):
        """Tests if self is compatible with other.

        Two CSEs are compatible if all they subexpressions are compatible.

        Args:
            other (CSE): The other CSE.

        Returns:
            True if self is compatible with other, False otherwise.
        """
        for self_subexpr, other_subexpr in itertools.product(self.subexprs, other.subexprs):
            if not self_subexpr.is_compatible(other_subexpr):
                return False
        return True

class CSEDetector():
    """A wrapper for a dictionary to detect common subexpressions.

    One of the challenges of detecting common subexpression is considering that
    for example A^T B = B^T A. The idea to solve this is to construct a
    surjective mapping of expressions to lists of Subexpression objects.

    Whenever a subexpression is added, for all variants of the subexpression
    (transposed, inverted) an entry is added to self.CSE_detection_dict that
    points to the same list of subexpressions, containing the newly added
    subexpression. Those lists of subexpressions are then potential common
    subexpressions.
    """
    def __init__(self):
        super().__init__()
        self.CSE_detection_dict = dict()
        self.all_subexpressions = dict()
        self.all_CSEs = dict()
        self.subexpr_to_CSE = dict()
        self.all_subexpr_lists = []

    def add_subexpression(self, subexpr):
        """Adds a subexpression to the detector.

        Args:
            subexpr (Subexpression): The subexpression to add.
        """
        for expr in subexpr.expr_var:
            try:
                subexprs = self.CSE_detection_dict[expr]
            except KeyError:
                subexprs = []
                self.all_subexpr_lists.append(subexprs)
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
        """Constructs CSE objects from all subexpressions that were added.

        The CSE objects are not returned, but stored in the self.all_CSEs
        dictionary.
        """
        for subexprs in self.all_subexpr_lists:
            if len(subexprs) > 1:
                cse = CSE(subexprs)
                # print(str(cse.expr), cse.max_clique_size)
                # Removing CSEs that can never be replaced.
                if cse.max_clique_size >= 2:
                    self.all_CSEs[cse.id] = cse
                    for subexpr in subexprs:
                        self.subexpr_to_CSE[subexpr.id] = cse.id

        # print("construct CSEs", [(str(cse.expr), len(cse.subexprs)) for cse in self.all_CSEs.values()])


    def maximal_CSEs(self):
        """Returns all maximal-replaceable CSEs.

        A common subexpression is maximal-replaceable if for all common
        subexpressions that subsume it, the number of replaceable occurrences
        is smaller.

        Returns:
            list: The list of maximal CSEs.
        """

        maximal_CSEs = []
        for cse in self.all_CSEs.values():
            if all(not cse.is_subexpression(_cse) or cse.max_clique_size > _cse.max_clique_size for _cse in self.all_CSEs.values()):
                maximal_CSEs.append(cse)

        return maximal_CSEs


    def replaceable_CSEs(self, CSEs):
        """Finds replaceable CSEs.
        
        A set of common subexpressions is replaceable if
        - every subexpression is compatible with every other subexpression
        - for each CSE, there are at least two subexpressions.

        Args:
            CSEs (list): A list of CSEs.

        Yields:
            list: Sets of replaceable CSEs, as a list of subexpressions.
        """
        if len(CSEs) == 1:
            # No constraints need to be checked here. All subexpressions belong
            # to the same CSE, and there are no cliques with less than two
            # nodes.
            for clique in CSEs[0].all_cliques:
                # TODO We yield maximal cliques here. Maximal cliques might
                # be smaller than maximum cliques. However, for the definition
                # of maximal CSEs, we rely on the clique number, that is, the
                # size of a maximum clique. So does it make sense to return
                # smaller cliques here? Probably not. However, maximal cliques
                # that are not maximum cliques are rare. The simplest case I can
                # come up with is X = AAAA, Y = AA.
                yield [self.all_subexpressions[id] for id in clique]
        else:
            subexprs = list(itertools.chain.from_iterable([node.subexprs for node in CSEs]))
            adj = construct_graph(subexprs)

            for clique in find_max_cliques(adj):
                # TODO Do maximal cliques that are not maximum cliques make
                # sense here?
                counter = Counter(self.subexpr_to_CSE[id] for id in clique)

                # Subsets are only used if there are at least two occurrences of
                # each CSE.
                if len(counter) == 2 and all(v >= 2 for v in counter.values()):
                    yield [self.all_subexpressions[id] for id in clique]

    def CSEs(self):
        """Generates all detected CSEs.

        This function generates the following replaceable common subexpressions:
        - Maximal common subexpressions.
        - Pairs of two maximal common subexpression where that are not
          subexpressions of each other and have overlapping occurrences.

        For all of those subexpressions, all subsets of subexpressions are
        generated which are replaceable.

        Yields:
            list: Sets of replaceable CSEs, as a list of subexpressions.
        """
        self.construct_CSEs()
        
        grouped_CSEs = []

        maximal_CSEs = self.maximal_CSEs()
        for cse in maximal_CSEs:
            grouped_CSEs.append((cse,))  

        # for cse1, cse2 in itertools.combinations(maximal_CSEs, 2):
        #     if (not (cse1.is_subexpression(cse2) or cse2.is_subexpression(cse1))) and not cse1.is_compatible(cse2):
        #         grouped_CSEs.append((cse1, cse2))

        for group in grouped_CSEs:
            # print("group", [str(cse.expr) for cse in group])
            yield from self.replaceable_CSEs(group)


def construct_graph(subexprs):
    """Constructs compatibility graph of subexpressions.

    Args:
        subexprs (list): List of subexpressions.

    Returns:
        dict: Maps subexpression IDs to sets of IDs of compatible
              subexpressions.
    """

    adj = dict()
    for subexpr1, subexpr2 in itertools.combinations(subexprs, 2):
        if subexpr1.is_compatible(subexpr2):
            adj.setdefault(subexpr1.id, set()).add(subexpr2.id)
            adj.setdefault(subexpr2.id, set()).add(subexpr1.id)
    return adj


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

def all_subexpressions(expr, position=tuple(), level=0, predeccessor=None):
    """Generates all subexpressions of expr.

    Generates all subexpression of an expression for CSE detection. Technically,
    not all subexpression are generated because it is already tested if
    expressions are blocked. Also, single operands, transposed operands, and
    expression that would introduce explicit inversion when replaced are not
    returned either.

    Canonical positions: If a subexpression is a single node in the expression
    tree, its position is only one sequence.

    Example: For ABC+D
    ABC+D   {()}                not {(0,), (1,)}
    ABC     {(0,)}              not {(0, 0), (0, 1), (0, 2)}
    D       {(1,)}
    AB      {(0, 0), (0, 1)}

    This way, positions are unique.

    Yields:
        Expression: The subexpression.
        set: A set of positions.
        int: The level of the expression. This is also the length of the positions.
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
                        new_expr = Plus(*_ops)
                        if not is_blocked(new_expr):
                            yield new_expr, positions, level+1
        elif expr.associative:
            if not is_blocked(expr):
                yield expr, {position}, level
                for i in range(2, len(expr.operands)):
                    for offset, seq in enumerate(window(expr.operands, i)):
                        positions = set(position + (offset+j,) for j in range(i))
                        new_expr = Times(*seq)
                        if not is_blocked(new_expr):
                            yield new_expr, positions, level+1
        for n, operand in enumerate(expr.operands):
            new_position = position + (n,)
            yield from all_subexpressions(operand, new_position, level+1, expr)

def indentify_subexpression_types(subexprs):
    """Assigns CSEType to each subexpression.

    For now, this function generates one arbitrary assignement.

    TODO what about A^-1 B^-1 ?

    Args:
        subexprs (list): A list of subexpressions.

    Returns:
        list: List of CSEType objects, one for each element of subexprs.
    """

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
    
    return return_types

def sort_keyfunc(CSE):
    """Keyfunction for sorting CSEs.

    As a sorting key for a common subexpression, this function returns a tuple
    of two elements, containing:
    - The overall number of symbols replaced, i.e. the number of occurrences
      times the number of symbols (operands) per occurrence.
    - The number of occurrences (to break ties).

    Returns:
        (int, int)
    """
    number_of_symbols = 0
    for expr, _ in CSE[0].expr.preorder_iter():
        if isinstance(expr, Symbol):
            number_of_symbols += 1

    return (len(CSE)*number_of_symbols, len(CSE))

def find_CSEs(equations):
    """Finds and replaces common subexpressions in equations.

    Args:
        equations (Equations): Some equations.

    Yields:
        Equations: The input equation with eliminated common subexpressions.        
    """

    CSE_detector = CSEDetector()
    for eqn_idx, equation in enumerate(equations):
        for expr, positions, level in all_subexpressions(equation.rhs):
            # print(expr, positions)

            # for expressions of the form Inverse(expr), we don't want to add the
            # inverted variant because this will produce "fake" CSEs. The reason
            # is that expr will also be added.
            # same for transpose
            # TODO do we want to invert if expr is not square?
            inv = contains_inverse(expr) and not is_inverse(expr)
            trans = contains_transpose(expr) and not is_transpose(expr)
            if expr.has_property(Property.SYMMETRIC):
                trans = False
            if inv and trans:
                subexpr = Subexpression(expr, eqn_idx, positions, level, CSEType.inverse_transpose)
                CSE_detector.add_subexpression(subexpr)
            elif inv:
                subexpr = Subexpression(expr, eqn_idx, positions, level, CSEType.inverse)
                CSE_detector.add_subexpression(subexpr)
            elif trans:
                subexpr = Subexpression(expr, eqn_idx, positions, level, CSEType.transpose)
                CSE_detector.add_subexpression(subexpr)
            else:
                subexpr = Subexpression(expr, eqn_idx, positions, level, CSEType.none)
                CSE_detector.add_subexpression(subexpr)

    CSEs = list(CSE_detector.CSEs())
    CSEs.sort(key=sort_keyfunc, reverse=True)

    for CSE in CSEs:
        # print("CSEs", [(str(subexpr.expr), subexpr.eqn_idx) for subexpr in sorted(CSE)])

        CSE_as_dict = dict()
        for subexpr in CSE:
            CSE_as_dict.setdefault(CSE_detector.subexpr_to_CSE[subexpr.id], []).append(subexpr)
  
        insert_equations = []
        replacements_per_equation = dict()
        for CSE_id, subexprs in CSE_as_dict.items():
            subexprs.sort()
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

        yield new_equations


if __name__ == "__main__":

    import linnea.config

    linnea.config.init()

    from linnea.derivation.graph.constructive import DerivationGraph
    # from linnea.derivation.graph.exhaustive import DerivationGraph
    # from linnea.derivation.graph.matrix_chain_derivation import DerivationGraph

    import linnea.examples.examples

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
