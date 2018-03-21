import string
import textwrap
import itertools

import matchpy


class CodeTemplate(object):
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
        if not isinstance(other, PropertyConstraint):
            return NotImplemented
        return self.variable == other.variable and self.properties == other.properties

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