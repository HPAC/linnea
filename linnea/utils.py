from .algebra import expression as ae
from .algebra import properties

import string
import textwrap
import itertools

import matchpy


def is_inverse(expr):
    return isinstance(expr, (ae.Inverse, ae.InverseTranspose, ae.InverseConjugate, ae.InverseConjugateTranspose))


def is_transpose(expr):
    return isinstance(expr, (ae.Transpose, ae.InverseTranspose, ae.ConjugateTranspose, ae.InverseConjugateTranspose))


def contains_inverse(expr):
    if is_inverse(expr):
        return True
    if isinstance(expr, ae.Operator):
        return any(contains_inverse(operand) for operand in expr.operands)


def contains_transpose(expr):
    if is_transpose(expr):
        return True
    if isinstance(expr, ae.Operator):
        return any(contains_transpose(operand) for operand in expr.operands)


class CodeTemplate():
    """Represents signatures of kernels.

    This class is mostly a wrapper for a string.Template object. The purpose is
    to allow for simpler repeated partial subsitutions and simple renaming of
    identifiers.

    Attributes:
        template (string.Template): Represents the signature.
    """
    def __init__(self, template=""):
        self.template = string.Template(textwrap.dedent(template))
    
    def __repr__(self):
        return self.template.safe_substitute()

    def safe_substitute(self, mapping=dict(), **kwds):
        """Safe substitutions in place.

        Note:
            This function modifies the object.
        """
        partial_substituted_str = self.template.safe_substitute(mapping, **kwds)
        self.template = string.Template(partial_substituted_str)

    def safe_substitute_copy(self, mapping=dict(), **kwds):
        """Safe substitutions returning a new CodeTemplate object.

        Note:
            This function does not modify the object.

        Returns:
            A copy of this CodeTemplate object with substitutions applied.
        """
        partial_substituted_str = self.template.safe_substitute(mapping, **kwds)
        return CodeTemplate(partial_substituted_str)

    def safe_substitute_str(self, mapping=dict(), **kwds):
        """Safe substitutions returning a string.

        Note:
            This function does not modify the object.

        Returns:
            A string with substitutions applied.
        """
        return self.template.safe_substitute(mapping, **kwds)

    def substitute_identifiers(self, mapping=dict(), **kwds):
        """Renames idendifiers in the signature.

        Args:
            mapping (dict): Dictionary that maps strings to strings. The keys
                are the old names of the identifiers and the values are the new
                names.
        """
        kwds.update(mapping)
        substitution = {key:"".join(["$", value]) for key, value in kwds.items()}
        self.safe_substitute(substitution)

    def substitute_identifiers_copy(self, mapping=dict(), **kwds):
        """Renames idendifiers in the signature, returning a copy.

        Args:
            mapping (dict): Dictionary that maps strings to strings. The keys
                are the old names of the identifiers and the values are the new
                names.
        """
        kwds.update(mapping)
        substitution = {key:"".join(["$", value]) for key, value in kwds.items()}
        return self.safe_substitute_copy(substitution)


def clamp(val, minval, maxval):
    """Returns minval for val < minval and maxval for val > maxval"""
    if val < minval: return minval
    if val > maxval: return maxval
    return val

def _sum(a, b):
    """Computes the sum of all integers from a to b.
    """
    return (b-a+1)*(a+b)//2

def _number_of_entries(bw, d1, d2):
    """Computes the number of nonzero entries contributed by bw.

    bw is the bandwidth. d1 and d2 are the matrix dimensions.

    bw has to be the bandwidth associated with the dimension d1, i.e.
    bw  d1   d2
    ------------
    lb rows cols
    ub cols rows
    """
    # For the comments below, assume that
    # bw = lb
    # d1 = n = rows
    # d2 = m = cols
    if d1 > d2:
        z = d1-d2
        # alpha is the number of full-length bands below the main diagonal
        # each full-length band has m entries
        alpha = min(bw, z)
        # beta is the number of non full-length bands
        # those bands have length m-1, m-2,…, m-beta
        beta = bw - alpha
        entries = alpha*d2 + _sum(d2-beta, d2-1)
    else:
        # there are no full-length bands in this case
        # bands have length n-1, n-2,…, n-bw
        entries = _sum(d1-bw, d1-1)
    return entries

def number_of_entries(expr):
    n, m = expr.size
    lb, ub = expr.bandwidth
    if lb >= 0 and ub >= 0:
        entries = min(n, m) # contribution of the diagonal
        entries += _number_of_entries(lb, n, m)
        entries += _number_of_entries(ub, m, n)
    elif lb < 0:
        entries = _number_of_entries(ub, m, n)
        entries -= _number_of_entries(-lb-1, m, n)
    elif ub < 0:
        entries = _number_of_entries(lb, n, m)
        entries -= _number_of_entries(-ub-1, n, m)
    if entries < 0:
        entries = 0
    return entries


######################
# partial ordering

def transitive_closure(closure):
    while True:
        new_relations = set((x,w) for x,y in closure for q,w in closure if q == y)
        closure_until_now = closure | new_relations
        if closure_until_now == closure:
            break
        closure = closure_until_now
    return closure

def transitive_reduction(closure):
    """
    https://en.wikipedia.org/wiki/Transitive_reduction
    """
    while True:
        new_relations = set((x,w) for x,y in closure for q,w in closure if q == y)
        closure_until_now = closure - new_relations
        if closure_until_now == closure:
            break
        closure = closure_until_now
    return closure

def to_dot_file(self):
    out = ["digraph G {", "ranksep=1;", "rankdir=TB;"]

    for element in self:
        out.append("{0} [shape=box];".format(element.name))

    for e1, e2 in self.__ordering__:
        out.append("{0} -> {1};".format(self(e2).name, self(e1).name))

    out.append("}")
    
    file_name = "partial_ordering.gv"
    output_file = open(file_name, "wt")
    output_file.write("\n".join(out)) 
    print("Output was saved in %s" % file_name)
    output_file.close()


def PartiallyOrderedEnum(enum):
    ordering = transitive_closure(enum.__ordering__)

    def lt(self, other):
        if self.__class__ is other.__class__:
            if (self.value, other.value) in ordering:
                return True
            return False
        return NotImplemented

    def le(self, other):
        if self.__class__ is other.__class__:
            if self == other:
                return True
            elif (self.value, other.value) in ordering:
                return True
            return False
        return NotImplemented

    enum.__lt__ = lt
    enum.__le__ = le
    enum.to_dot_file = to_dot_file
    enum.__transitive_closure__ = ordering
    return enum


class PropertyConstraint(matchpy.Constraint):
    def __init__(self, variable, properties):
        self.variable = variable # variable name (str)
        # Removing redundant properties.
        remove = set()
        for p1, p2 in itertools.combinations(properties, 2):
            if p1 < p2:
                remove.add(p2)
                continue
            if p2 < p1:
                remove.add(p1)
        self.properties = frozenset(properties.difference(remove))

    def __call__(self, match):
        try:
            variable = match[self.variable]
        except KeyError:
            return False
        else:
            return all(variable.has_property(prop) for prop in self.properties)

    def __eq__(self, other):
        return isinstance(other, PropertyConstraint) and self.variable == other.variable and self.properties == other.properties

    def __hash__(self):
        return hash((self.variable, self.properties))

    def __repr__(self):
        return 'PropertyConstraint({!r}, {{{}}})'.format(self.variable, ', '.join(repr(p) for p in self.properties))

    def __str__(self):
        return 'PC({} is {})'.format(self.variable, ', '.join(str(p) for p in self.properties))

    def with_renamed_vars(self, renaming):
        try:
            return PropertyConstraint(renaming[self.variable], self.properties)
        except KeyError:
            return self

    @property
    def variables(self):
        return frozenset([self.variable])


class InequalityConstraint(matchpy.Constraint):
    def __init__(self, x1, x2):
        self.x1 = x1 # variable name (str)
        self.x2 = x2 # variable name (str)

    def __call__(self, match):
        try:
            x1 = match[self.x1]
        except KeyError:
            return False
        else:
            try:
                x2 = match[self.x2]
            except KeyError:
                return False
            else:
                return x1 != x2

    def __eq__(self, other):
        return isinstance(other, InequalityConstraint) and self.x1 == other.x1 and self.x2 == other.x2

    def __hash__(self):
        return hash((self.x1, self.x2))
 
    def __repr__(self):
        return 'InequalityConstraint({!r}, {!r})'.format(self.x1, self.x2)

    def __str__(self):
        return 'IC({} != {})'.format(self.x1, self.x2)

    def with_renamed_vars(self, renaming):
        try:
            return InequalityConstraint(renaming[self.x1], renaming[self.x2])
        except KeyError:
            return self

    @property
    def variables(self):
        return frozenset([self.x1, self.x2])


class Element(object):
    """Element for UnionFind."""
    def __init__(self, elem):
        super(Element, self).__init__()
        self.content = elem
        self.parent = self
        self.rank = 0
        

class UnionFind(object):
    """docstring for UnionFind"""
    def __init__(self):
        super(UnionFind, self).__init__()
        self.elements = dict()
        self.components = dict()

    def add_element(self, elem):
        if not elem in self.elements:
            obj = Element(elem)
            # self.content.append(obj)
            self.elements[elem] = obj
            self.components[elem] = [obj]

    # def find(self, elem):
    #     root = self.elements[elem]
    #     while root.parent != root:
    #         root = root.parent

    #     # while elem.parent != root:
    #     #     parent = elem.parent
    #     #     elem.parent = root
    #     #     elem = parent

    #     return root

    def find(self, elem):
        """

        Implements the path halving optimization.
        """
        node = self.elements[elem]
        while node.parent != node:
            node.parent = node.parent.parent
            node = node.parent
        return node
        
    def union(self, elem1, elem2):
        root1 = self.find(elem1)
        root2 = self.find(elem2)

        if root1 == root2:
            return

        if root1.rank < root2.rank:
            root1, root2 = root2, root1

        self.components[root1.content].extend(self.components[root2.content])
        del self.components[root2.content]

        root2.parent = root1
        if root1.rank == root2.rank:
            root1.rank = root1.rank + 1

    def sets(self):
        return [set(elem.content for elem in components) for components in self.components.values()]
            

def collect_dependent_dimensions(expr, union_find):
    if isinstance(expr, ae.Symbol):
        rows = (expr.name, 0)
        cols = (expr.name, 1)
        union_find.add_element(rows)
        union_find.add_element(cols)
        if expr.has_property(properties.Property.SYMMETRIC) and not expr.has_property(properties.Property.DIAGONAL):
            union_find.union(rows, cols)
        return (rows, cols)
    elif isinstance(expr, ae.Transpose):
        rows, cols = collect_dependent_dimensions(expr.operand, union_find)
        return (cols, rows)
    elif isinstance(expr, ae.InverseTranspose):
        rows, cols = collect_dependent_dimensions(expr.operand, union_find)
        union_find.union(rows, cols) # inverted expressions have to be square
        return (rows, cols)
    elif isinstance(expr, ae.Inverse):
        rows, cols = collect_dependent_dimensions(expr.operand, union_find)
        union_find.union(rows, cols) # inverted expressions have to be square
        return (cols, rows)
    elif isinstance(expr, ae.Times):
        dimension_list = [collect_dependent_dimensions(op, union_find) for op in expr.operands]
        stack = []
        n = len(expr.operands)
        rows = None
        cols = None
        for i in range(1, n):
            current_op = expr.operands[i]
            previous_op = expr.operands[i-1]
            if current_op.rows != 1:
                if previous_op.columns != 1:
                    # no scalars
                    union_find.union(dimension_list[i-1][1], dimension_list[i][0])
                else:
                    # end of scalar subexpression
                    if stack:
                        last_dim = stack.pop()
                        union_find.union(last_dim, dimension_list[i][0])
                    else:
                        # expressions starts with scalar subexpression
                        rows = dimension_list[i][0]
            else:
                if previous_op.columns != 1:
                    # beginning of scalar subexpression
                    stack.append(dimension_list[i-1][1])
                # else: consecutive scalar subexpressions

        if stack:
            # expression ends with scalar subexpression
            cols = stack.pop()
        else:
            # expression does not end with scalar subexpression
            cols = dimension_list[-1][1]
        if not rows:
            # expression does not start with scalar subexpression
            rows = dimension_list[0][0]

        return (rows, cols)
    elif isinstance(expr, ae.Plus):
        dimension_list = [collect_dependent_dimensions(op, union_find) for op in expr.operands]
        for op1, op2 in window(dimension_list):
            union_find.union(op1[0], op2[0])
            union_find.union(op1[1], op2[1])
        return dimension_list[0]
    elif isinstance(expr, ae.Equal):
        rhs_dimensions = collect_dependent_dimensions(expr.rhs, union_find)
        lhs_dimensions = collect_dependent_dimensions(expr.lhs, union_find)
        union_find.union(rhs_dimensions[0], lhs_dimensions[0])
        union_find.union(rhs_dimensions[1], lhs_dimensions[1])
        return rhs_dimensions


def dependent_dimensions(equations):
    """Computes dependent dimensions.

    The dependent dimensions are all sets of dimensions that have to be the same
    for the input equations to be valid. Dimensions are represented as tuples of
    two elements: The first element is the name of the operand, the second
    element is an integer; 0 stands for rows, 1 for columns.

    Args:
        equations (Equations): Input equations

    Returns:
        list: A list of sets. All dimensions in one set have to be the same.
    """
    uf = UnionFind()
    for eqn in equations:
        collect_dependent_dimensions(eqn, uf)
    return uf.sets()


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


def roundrobin(*iterables):
    "roundrobin('ABC', 'D', 'EF') --> A D E B F C"
    # Recipe credited to George Sakkis
    num_active = len(iterables)
    nexts = itertools.cycle(iterables)
    while num_active:
        try:
            for n in nexts:
                yield next(n)
        except StopIteration:
            # Remove the iterator we just exhausted from the cycle.
            num_active -= 1
            nexts = itertools.cycle(itertools.islice(nexts, num_active))