import matchpy
import functools
import operator

from . import InferenceOfProperties
from .properties import Property as properties
from .properties import implications
from . import utils

from .. import config

from ..config import CppLibrary

dot_table_node = """{0} [shape=record, label="{{{{ {1} | {2} }}| {3} | {4} }}"];\n"""

class Expression(matchpy.Expression):
    """docstring for Expression"""

    known_inv_symbols = set() # necessary for computing metric
    known_symbols = set() # necessary for computing metric

    counter = 0 # for to_dot()

    def __init__(self, variable_name=None):
        super(Expression, self).__init__(variable_name)

    def init_partitioning(self):
        """Adds the attribute "partitioning" to all nodes in the expression tree.

        Copying the partitioning attribute (a tuple of two sets) is surprisingly
        expensive. However, this attribute is only needed once when the
        partitioning is distributed. Thus, to speed up copying expression,
        this attribute can be set and removed during runtime. The custom
        implementation of deepcopy does not copy this attribute.
        """
        # partitioning is only set if it doesn't exist yet to avoid overwriting
        # partitionings defined by the user
        try:
            self.partitioning
        except AttributeError:
            self.partitioning = (set(), set())

        if not isinstance(self, BlockedExpression) and not isinstance(self, Symbol):
            for operand in self.operands:
                operand.init_partitioning()

    def delete_partitioning(self):
        """Removes the attribute "partitioning" to all nodes in the expression tree.

        See docstring for init_partitioning().
        """
        # This is necessary because it's possible that symbols are used multiple
        # times in multiple expressions. Thus, this function may visit a symbol
        # that was already visited before.
        try:
            del self.partitioning
        except AttributeError:
            pass
        # TODO do I need to check for BlockeExpressions here, too? This function
        # should probably be called when there are not BlockedExpressions anymore
        if not isinstance(self, Symbol):
            for operand in self.operands:
                operand.delete_partitioning()

    def has_property(self, prop):
        return InferenceOfProperties.infer_property(self, prop)

    def index_range(self):
        return functools.reduce(operator.mul, [index.value for index in self.indices], 1)

    def metric(self):
        """Returns the metric of self

        The metric is a list of two integers. The first one is a count of the
        number of unique (i.e. duplicates are not counted) matrices that may
        be considered for factoring (matrices that are not triangular, diagonal
        or orthogonal). The second integer is a count of all nodes (i.e.
        operands) in the expression (counting duplicates TODO change that?).
        """
        Expression.known_inv_symbols = set()
        Expression.known_symbols = set()
        return self._metric()
        # m = self._metric_unique()
        # print(self)
        # print(Expression.known_inv_symbols)
        # print(Expression.known_symbols)
        # return m
        # return self._metric_unique()

    def _metric(self, inverse=False):
        metric = [0, 0]
        if isinstance(self, Symbol):
            metric[1] += 1
            if inverse and self.name not in Expression.known_inv_symbols:
                Expression.known_inv_symbols.add(self.name)
                if not any(self.has_property(prop) for prop in [properties.TRIANGULAR, properties.DIAGONAL, properties.ORTHOGONAL]):
                    metric[0] += 1
        else:
            inv = (isinstance(self, Inverse) or isinstance(self, InverseTranspose) or isinstance(self, InverseConjugate) or isinstance(self, InverseConjugateTranspose) or inverse)
            for operand in self.operands:
                operand_metric = operand._metric(inv)
                metric[0] += operand_metric[0]
                metric[1] += operand_metric[1]
        return metric


    def _metric_unique(self, inverse=False):
        # This funciton just calls unique nodes.
        metric = [0, 0]
        if isinstance(self, Symbol):
            # metric[1] += 1
            if inverse and self.name not in Expression.known_inv_symbols:
                Expression.known_inv_symbols.add(self.name)
                if not any(self.has_property(prop) for prop in [properties.TRIANGULAR, properties.DIAGONAL, properties.ORTHOGONAL]):
                    metric[0] += 1
            if self.name not in Expression.known_symbols:
                Expression.known_symbols.add(self.name)
                metric[1] += 1
        else:
            inv = (isinstance(self, Inverse) or isinstance(self, InverseTranspose) or isinstance(self, InverseConjugate) or isinstance(self, InverseConjugateTranspose) or inverse)
            for operand in self.operands:
                operand_metric = operand._metric(inv)
                metric[0] += operand_metric[0]
                metric[1] += operand_metric[1]
        return metric

        # self_is_inv = isinstance(self, Inverse)

        # metric = 0
        # for ch in self.operands:
        #     if isinstance(ch, Symbol):
        #         # alternatively all(not has_property)
        #         if not any([ch.has_property(prop) for prop in [properties.TRIANGULAR, properties.DIAGONAL, properties.ORTHOGONAL]]):
        #             metric += 1
        #     else:
        #         metric += ch.metric(self_is_inv)
        # return metric


    def inverse_of(self, other):
        """Tests if self == other^-1 holds.

        For this to work, both self and other have to be simplified.
        """
        if isinstance(self, Inverse):
            return self.operand == other
        elif isinstance(other, Inverse):
            return self == other.operand
        elif isinstance(self, InverseTranspose) and isinstance(other, Transpose):
            return self.operand == other.operand
        elif isinstance(self, Transpose) and isinstance(other, InverseTranspose):
            return self.operand == other.operand
        elif isinstance(self, InverseConjugate) and isinstance(other, Conjugate):
            return self.operand == other.operand
        elif isinstance(self, Conjugate) and isinstance(other, InverseConjugate):
            return self.operand == other.operand
        elif isinstance(self, InverseConjugateTranspose) and isinstance(other, ConjugateTranspose):
            return self.operand == other.operand
        elif isinstance(self, ConjugateTranspose) and isinstance(other, InverseConjugateTranspose):
            return self.operand == other.operand
        elif isinstance(self, Symbol) and self.has_property(properties.ORTHOGONAL) and isinstance(other, Transpose) and other.operand.has_property(properties.ORTHOGONAL):
            return self == other.operand
        elif isinstance(other, Symbol) and other.has_property(properties.ORTHOGONAL) and isinstance(self, Transpose) and self.operand.has_property(properties.ORTHOGONAL):
            return other == self.operand
        elif isinstance(self, Symbol) and isinstance(other, Symbol):
            if self.has_property(properties.IDENTITY) and self == other:
                return True
            else:
                return False
        # the inverse_of relation is not affected by Transpose, ConjugateTranspose and Conjugate
        # TODO use set instead of list? Construct this list/set just once?
        elif type(self) == type(other) and type(self) in [Transpose, ConjugateTranspose, Conjugate]:
            return self.operand.inverse_of(other.operand)
        elif isinstance(self, Times) and isinstance(other, Times):
            # print("HERE", self, other)
            self_scalars, self_non_scalars = self.split_operands()
            other_scalars, other_non_scalars = other.split_operands()

            # print(self_scalars, other_scalars)
            # print(self_non_scalars, other_non_scalars)

            inv_of = True
            for s, o in zip(_inverse_iterator_LR(self_non_scalars), _inverse_iterator_RL(other_non_scalars)):
                self_operand, self_inside_inverse = s
                other_operand, other_inside_inverse = o
                if self_inside_inverse or other_inside_inverse:
                    inv_of = (self_operand == other_operand)
                else:
                    inv_of = self_operand.inverse_of(other_operand)
                if not inv_of:
                    return False

            return _test_relation_commutative(self_scalars, other_scalars, "inverse_of")
        elif isinstance(self, Plus) and isinstance(other, Plus):
            return _test_relation_commutative(self.operands, other.operands, "inverse_of")
        else:
            return False

    def transpose_of(self, other):
        """Returns true if self is the transpose of other.

        To ensure that this function returns a correct results, is is necessary
        that both expressions are in the same canonical form. It should always
        return a correct result if self == other (i.e. when testing symmetry).
        """
        if self.has_property(properties.SCALAR) and other.has_property(properties.SCALAR):
            return self == other
        elif isinstance(self, Transpose):
            return self.operand == other
        elif isinstance(other, Transpose):
            return self == other.operand
        elif isinstance(self, InverseTranspose) and isinstance(other, Inverse):
            return self.operand == other.operand
        elif isinstance(self, Inverse) and isinstance(other, InverseTranspose):
            return self.operand == other.operand
        elif isinstance(self, ConjugateTranspose) and isinstance(other, Conjugate):
            return self.operand == other.operand
        elif isinstance(self, Conjugate) and isinstance(other, ConjugateTranspose):
            return self.operand == other.operand
        elif isinstance(self, InverseConjugateTranspose) and isinstance(other, InverseConjugate):
            return self.operand == other.operand
        elif isinstance(self, InverseConjugate) and isinstance(other, InverseConjugateTranspose):
            return self.operand == other.operand
        elif isinstance(self, Symbol) and isinstance(other, Symbol):
            if self.has_property(properties.SYMMETRIC) and self == other:
                return True
            else:
                return False
        # the transpose_of relation is not affected by Inverse, InverseConjugate and Conjugate
        # TODO use set instead of list? Construct this list/set just once?
        elif type(self) == type(other) and type(self) in [Inverse, InverseConjugate, Conjugate]:
            return self.operand.transpose_of(other.operand)
        elif isinstance(self, Times) and isinstance(other, Times):
            # print("times", list(zip(self.operands, reversed(other.operands))))
            if len(self.operands) != len(other.operands):
                return False
            self_scalars, self_non_scalars = self.split_operands()
            other_scalars, other_non_scalars = other.split_operands()
            return all(self_operand.transpose_of(other_operand) for self_operand, other_operand in zip(self_non_scalars, reversed(other_non_scalars))) and _test_relation_commutative(self_scalars, other_scalars, "transpose_of")
        elif isinstance(self, Plus) and isinstance(other, Plus):
            return _test_relation_commutative(self.operands, other.operands, "transpose_of")
        else:
            return False

    def to_dot_file(self, file_name=None):

        Expression.counter = 0

        out = ""
        out = self.to_dot(out)
        out = "\n".join(["digraph G {", "ranksep=1.5;", "ordering=out;", "rankdir=TB;", out, "}"])

        if not file_name:
            file_name = str(self) + ".gv"
        file_name = "debug_output/expressions/" + file_name
        output_file = open(file_name, "w")
        output_file.write(out)
        output_file.close()
        print("Output was saved in %s" % file_name)

    def to_julia_expression(self):
        raise NotImplementedError()

def _test_relation_commutative(l1, l2, relation):
    """Helper function for inverse_of, transpose_of

    When testing if Plus(A, B^T) is the transpose of Plus(B^T, A), the fact that
    Plus is commutative causes a problem: The positions of summands in one
    expression and transposed summands in the other expression are not known.
    Sorting the summands doesn't help because they would be sorted
    differently. So it is necessary to search for the corresponding transposed
    terms.

    This is what this function does. It takes two list as input. The third
    argument is the name of the relation that is tested. This function returns
    True if for each element in l1, there is one element in l2 sucht that the
    relation holds for those two elements.

    The function with the name given in the argument relation is called as a
    method on elements of l1 (i.e. l1[0].relation(l2[0])), so this function has
    to be a method of all those elements.

    This is relevant for transpose_of and inverse_of. This problem also occurs
    with scalar products.
    """
    checked = set()
    for element1 in l1:
        for n, element2 in enumerate(l2):
            if n not in checked:
                if getattr(element1, relation)(element2):
                # if self_child.relation(element2):
                    checked.add(n)
                    break
        else:
            # When no corresponding element2 for a given element1 is
            # found, the relation can not hold
            return False
    return True

class Operator(matchpy.Operation, Expression):
    """docstring for Operator"""
    def __init__(self, *operands, variable_name=None):
        super(Operator, self).__init__(*operands, variable_name=variable_name)

    def split_operands(self):
        """Splits operands into scalar and non-scalar operands.

        Returns a pair of lists:
        - The first list contains all scalar operands.
        - The second list contains all non-scalar operands.

        This function does no change the order of the operands.
        """
        scalars = []
        non_scalars = []
        for operand in self.operands:
            try:
                is_scalar = operand.has_property(properties.SCALAR)
            except AttributeError:
                # This is the case if operand is a matchpy.Variable
                is_scalar = False
            if is_scalar:
                scalars.append(operand)
            else:
                non_scalars.append(operand)
        return scalars, non_scalars

    @property
    def indices(self):
        # this doesn't make sense for Equal
        return set.union(*[op.indices for op in self.operands])

    def to_dot(self, out):
        node_name = "".join(["node", str(Expression.counter)])
        properties = ""#"\n".join([str(p) for p in self.properties])
        false_properties = ""#"\n".join([str(p) for p in self.false_properties])
        out += dot_table_node.format(node_name, self.__class__.__name__, str(self.size), properties, false_properties)
        # out = "".join([out, node_name, " [shape=record, label=\"<f0>", str(self.__class__.__name__), "|<f1>", str(self.size), "|<f2>", str(self.properties), "|<f3>", str(self.false_properties), "\"];\n" ])
        # output_file.write(out)
        for operand in self.operands:
            Expression.counter+=1
            operand_name = "".join(["node", str(Expression.counter)])
            out = "".join([out, node_name, " -> ", operand_name, ";\n"])
            out = operand.to_dot(out)
        return out


class Symbol(matchpy.Symbol, Expression):
    """docstring for Symbol"""
    def __init__(self, name, size, indices):
        super(Symbol, self).__init__(name, variable_name=None)
        self.size = size
        self.indices = indices
        self.bandwidth = (size[0]-1, size[1]-1)
        self.properties = set()
        self.false_properties = set()

    def set_property(self, prop):
        # print("Setting property ", prop, " to ", self)

        if prop == properties.LOWER_TRIANGULAR:
            self.bandwidth = (self.bandwidth[0], 0)
        elif prop == properties.UPPER_TRIANGULAR:
            self.bandwidth = (0, self.bandwidth[1])
        elif prop == properties.DIAGONAL:
            self.bandwidth = (0, 0)

        # TODO those properties should not added to the sets because that causes
        # problems when using update_properties()
        if prop == properties.SCALAR:
            return
        elif prop == properties.VECTOR:
            return
        elif prop == properties.MATRIX:
            return

        self.properties.add(prop)
        self.properties.update(implications.get(prop, tuple()))

        self.false_properties.discard(prop)
        self.false_properties.difference_update(implications.get(prop, tuple()))

    def to_dot(self, out):
        node_name = "".join(["node", str(Expression.counter)])
        props = "\n".join([str(p) for p in self.properties])
        props = "\n".join([str(p) for p in properties if self.has_property(p)])
        false_properties = "\n".join([str(p) for p in self.false_properties])
        out += dot_table_node.format(node_name, str(self.name), str(self.size), props, false_properties)
        # out = "".join([out, node_name, " [shape=record, label=\"<f0>", str(self.name), "|<f1>", str(self.size), "|<f2>", str(self.properties), "|<f3>", str(self.false_properties), "\"];\n" ])
        # output_file.write(out)
        return out

    def to_julia_expression(self, recommended=False):
        return self.name

    def to_cpp_expression(self, lib, recommended=False):
        return self.name

class Constant(object):
    """docstring for Constant"""
    pass


######################
# operators
#


class Equal(Operator):
    """docstring for Equal"""

    name = "="
    arity = matchpy.Arity.binary
    associative = False
    commutative = False
    one_identity = False
    infix = True

    def __init__(self,  *operands, variable_name=None):
        super(Equal, self).__init__(*operands, variable_name=variable_name)

    @property
    def lhs(self):
        return self.operands[0]

    # @lhs.setter
    # def lhs(self, value):
    #     self.operands[0] = value

    @property
    def rhs(self):
        return self.operands[1]

    # @rhs.setter
    # def rhs(self, value):
    #     self.operands[1] = value

    @property
    def size(self):
        return self.operands[0].size

    @property
    def rows(self):
        return self.operands[0].rows

    @property
    def columns(self):
        return self.operands[0].columns

    def __str__(self):
        return "{0} = {1}".format(*self.operands)

    def to_julia_expression(self, recommended=False):
        return "{0} = {1}".format(*map(operator.methodcaller("to_julia_expression", recommended), self.operands))

    def to_cpp_expression(self, lib, recommended=False):
        return {
            CppLibrary.Blaze : "auto {0} = blaze::eval({1});",
            CppLibrary.Eigen : "auto {0} = ({1}).eval();",
            CppLibrary.Armadillo : "auto {0} = ({1}).eval();"
        }.get(lib).format(*map(operator.methodcaller("to_cpp_expression", lib, recommended), self.operands))

class Plus(Operator):
    """docstring for Plus"""

    name = "+"
    arity = matchpy.Arity.variadic
    associative = True
    commutative = True
    one_identity = True
    infix = True

    def __init__(self,  *operands, variable_name=None):
        super(Plus, self).__init__(*operands, variable_name=variable_name)

    @property
    def size(self):
        return self.operands[0].size

    @property
    def rows(self):
        return self.operands[0].rows

    @property
    def columns(self):
        return self.operands[0].columns

    @property
    def bandwidth(self):
        bands = list(zip(*[operand.bandwidth for operand in self.operands]))
        return (max(bands[0]), max(bands[1]))

    def to_julia_expression(self):
        operand_str = '+'.join(map(operator.methodcaller("to_julia_expression"), self.operands))
        return "({0})".format(operand_str)


class Times(Operator):
    """docstring for Times"""

    name = ""
    arity = matchpy.Arity.variadic
    associative = True
    commutative = False
    one_identity = True
    infix = True

    def __init__(self,  *operands, variable_name=None):
        super(Times, self).__init__(*operands, variable_name=variable_name)

    @property
    def size(self):
        # return (self.rows, self.columns)
        return self._calc_size()

    # @profile
    def _calc_size(self):
        scalar_expressions = utils.scalar_subexpressions(self)

        rows = 0
        cols = 0
        factors = self.operands
        length = len(factors)
        for pair in scalar_expressions:
            if pair[0] == 0:
                pos = pair[1]
                if pos == length-1:
                    rows = 1
                else:
                    rows = factors[pos+1].size[0]
            if pair[1] == length-1:
                pos = pair[0]
                if pos == 0:
                    cols = 1
                else:
                    cols = factors[pos-1].size[1]

        if rows == 0:
            rows = factors[0].size[0]
        if cols == 0:
            cols = factors[-1].size[1]

        return (rows, cols)


    @property
    # @profile
    def rows(self):
        stack = []
        factors = self.operands
        pos = -1 # With this initialization, even products with one operand work.
        for n, factor in enumerate(factors):
            if factor.rows == 1:
                stack.append(n)
            if factor.columns == 1 and stack:
                stack.pop()
                pos = n
            if not stack:
                break

        if pos == len(factors)-1:
            return 1
        else:
            return factors[pos+1].rows

    @property
    def columns(self):
        stack = []
        factors = self.operands
        pos = len(factors) # With this initialization, even products with one operand work.
        for n, factor in zip(itertools.count(pos-1, -1), reversed(self.operands)):
            if factor.columns == 1:
                stack.append(n)
            if factor.rows == 1 and stack:
                stack.pop()
                pos = n
            if not stack:
                break

        if pos == 0:
            return 1
        else:
            return factors[pos-1].columns

    @property
    def bandwidth(self):

        # Even though scalars have bandwidth (0, 0), which does not affect the
        # result, their size (1, 1) can affect the result. Thus, it's necessary
        # to use split_operands here.
        scalars, non_scalars = self.split_operands()
        if non_scalars:
            sizes = [operand.size for operand in non_scalars]
            bands = [operand.bandwidth for operand in non_scalars]
            l = len(non_scalars)
            # print(self, bands, sizes)
            lb = bands[-1][0]
            for i in range(l-2, -1, -1):
                lb = min(lb + bands[i][0], sizes[i][0]-1)
                if -sizes[i+1][1] >= lb:
                    lb = -self.size[1]
                    break

            ub = bands[0][1]
            for i in range(1, l):
                ub = min(ub + bands[i][1], sizes[i][1]-1)
                if -sizes[i-1][0] >= ub:
                    ub = -self.size[0]
                    break
            # print((lb, ub))
            return (lb, ub)
        else:
            return (0, 0)

    def __str__(self):
        operand_str = ' '.join(map(str, self.operands))
        return "({0})".format(operand_str)

    def is_inverse_type(self, type):
        return type in [Inverse, InverseTranspose]

    def recommended_julia_expression(self, idx):
        op = self.operands[idx]
        if self.is_inverse_type(type(op)):
            # there are at least two elements to process
            #  since op is inverse, it will be inv(a) * b -> a \ b
            # doesn't matter if b is inverse or not
            if idx != len(self.operands) - 1:
                return "({0}\{1})".format(
                    op.to_julia_expression(recommended=True, strip_inverse=True),
                    self.recommended_julia_expression(idx + 1)
                )
            else:
                return ("{0}*" if idx != len(self.operands) - 1 else "{0}").format(op.to_julia_expression(recommended=True))
        else:
            # only two elements are left
            if idx == len(self.operands) - 2:
                op_next = self.operands[idx + 1]
                # if next is an inverse, match a * inv(b) -> a / b
                if self.is_inverse_type(type(op_next)):
                    return "({0}/{1})".format(
                        op.to_julia_expression(recommended=True),
                        op_next.to_julia_expression(recommended=True, strip_inverse=True)
                    )
                # otherwise, simply finish
                else:
                    return "{0}*{1}".format(
                        op.to_julia_expression(recommended=True),
                        op_next.to_julia_expression(recommended=True)
                    )
            elif idx == len(self.operands) - 1:
                return op.to_julia_expression(recommended=True)
            else:
                return "{0}*{1}".format(
                    op.to_julia_expression(recommended=True),
                    self.recommended_julia_expression(idx + 1)
                )

    def recommended_cpp_expression(self, lib, idx):
        op = self.operands[idx]
        if self.is_inverse_type(type(op)):
            # there are at least two elements to process
            #  since op is inverse, it will be inv(a) * b -> a \ b
            # doesn't matter if b is inverse or not
            if idx != len(self.operands) - 1:
                return {
                    CppLibrary.Eigen: "({0}.partialPivLU().solve({1}))",
                    CppLibrary.Armadillo: "arma::solve({0}, {1}, arma::solve_opts::fast)"
                }.get(lib).format(
                    op.to_cpp_expression(lib, recommended=True, strip_inverse=True),
                    self.recommended_cpp_expression(lib, idx + 1)
                )
            else:
                return ("{0}*" if idx != len(self.operands) - 1 else "{0}").format(
                    op.to_cpp_expression(lib, recommended=True))
        else:
            # only two elements are left
            # leave if we decide to perform a/b for cpp as well
            if idx == len(self.operands) - 2:
                op_next = self.operands[idx + 1]
                return "{0}*{1}".format(
                    op.to_cpp_expression(lib, recommended=True),
                    op_next.to_cpp_expression(lib, recommended=True)
                )
            elif idx == len(self.operands) - 1:
                return op.to_cpp_expression(lib, recommended=True)
            else:
                return "{0}*{1}".format(
                    op.to_cpp_expression(lib, recommended=True),
                    self.recommended_cpp_expression(lib, idx + 1)
                )

    def to_julia_expression(self, recommended=False):
        # TODO are parenthesis necessary?
        # operand_str = '*'.join(map(operator.methodcaller("to_julia_expression"), self.operands))
        # return "({0})".format(operand_str)
        if recommended:
            return self.recommended_julia_expression(0)
        else:
            return '*'.join(map(operator.methodcaller("to_julia_expression"), self.operands))

    def to_cpp_expression(self, lib, recommended=False):
        if recommended:
            return self.recommended_cpp_expression(lib, 0)
        else:
            return '*'.join(map(operator.methodcaller("to_cpp_expression", lib), self.operands))


class LinSolveL(Operator):
    """docstring for LinSolveL

    This operator is only needed to generate Julia code.
    """

    name = "\\"
    arity = matchpy.Arity.binary
    associative = False
    commutative = False
    one_identity = False
    infix = True

    def to_julia_expression(self):
        return "({0}\{1})".format(*map(operator.methodcaller("to_julia_expression"), self.operands))


class LinSolveR(Operator):
    """docstring for LinSolveR

    This operator is only needed to generate Julia code.
    """

    name = "/"
    arity = matchpy.Arity.binary
    associative = False
    commutative = False
    one_identity = False
    infix = True

    def to_julia_expression(self):
        return "({0}/{1})".format(*map(operator.methodcaller("to_julia_expression"), self.operands))


class Identity(Operator):
    """docstring for Identity"""

    name = "Id"
    arity = matchpy.Arity.unary
    associative = False
    commutative = False
    one_identity = False
    infix = False

    def __init__(self,  *operands, variable_name=None):
        super(Identity, self).__init__(*operands, variable_name=variable_name)

    @property
    def size(self):
        return self.operands[0].size

    @property
    def rows(self):
        return self.operands[0].rows

    @property
    def columns(self):
        return self.operands[0].columns

    @property
    def operand(self):
        return self.operands[0]

    @operand.setter
    def operand(self, value):
        self.operands[0] = value

    @property
    def bandwidth(self):
        return self.operands[0].bandwidth

    def __str__(self):
        return "{1}({0})".format(self.operands[0], self.name)

class Transpose(Operator):
    """docstring for Transpose"""

    name = "^T"
    arity = matchpy.Arity.unary
    infix = False

    def __init__(self,  *operands, variable_name=None):
        super(Transpose, self).__init__(*operands, variable_name=variable_name)

    @property
    def size(self):
        return tuple(reversed(self.operands[0].size))

    @property
    def rows(self):
        return self.operands[0].columns

    @property
    def columns(self):
        return self.operands[0].rows

    @property
    def operand(self):
        return self.operands[0]

    @operand.setter
    def operand(self, value):
        self.operands[0] = value

    @property
    def bandwidth(self):
        return tuple(reversed(self.operands[0].bandwidth))

    def __str__(self):
        return "{0}{1}".format(self.operands[0], self.name)

    def to_julia_expression(self, recommended=False):
        return "{0}'".format(self.operands[0].to_julia_expression(recommended))

    def to_cpp_expression(self, lib, recommended=False):
        op = self.operands[0].to_cpp_expression(lib, recommended)
        return {
            CppLibrary.Blaze: "blaze::trans({0})",
            CppLibrary.Eigen: "({0}).transpose()",
            CppLibrary.Armadillo: "({0}).t()"
        }.get(lib).format(op)


class Conjugate(Operator):
    """docstring for Conjugate"""

    name = "^C"
    arity = matchpy.Arity.unary
    infix = False

    def __init__(self,  *operands, variable_name=None):
        super(Conjugate, self).__init__(*operands, variable_name=variable_name)

    @property
    def size(self):
        return self.operands[0].size

    @property
    def rows(self):
        return self.operands[0].rows

    @property
    def columns(self):
        return self.operands[0].columns

    @property
    def operand(self):
        return self.operands[0]

    @operand.setter
    def operand(self, value):
        self.operands[0] = value

    @property
    def bandwidth(self):
        return self.operands[0].bandwidth

    def __str__(self):
        return "{0}{1}".format(self.operands[0], self.name)


class ConjugateTranspose(Operator):
    """docstring for ConjugateTranspose"""

    name = "^H"
    arity = matchpy.Arity.unary
    infix = False

    def __init__(self,  *operands, variable_name=None):
        super(ConjugateTranspose, self).__init__(*operands, variable_name=variable_name)

    @property
    def size(self):
        return tuple(reversed(self.operands[0].size))

    @property
    def rows(self):
        return self.operands[0].columns

    @property
    def columns(self):
        return self.operands[0].rows

    @property
    def operand(self):
        return self.operands[0]

    @operand.setter
    def operand(self, value):
        self.operands[0] = value

    @property
    def bandwidth(self):
        return tuple(reversed(self.operands[0].bandwidth))

    def __str__(self):
        return "{0}{1}".format(self.operands[0], self.name)


class Inverse(Operator):
    """docstring for Inverse"""

    name = "^-1"
    arity = matchpy.Arity.unary
    infix = False

    def __init__(self,  *operands, variable_name=None):
        super(Inverse, self).__init__(*operands, variable_name=variable_name)

    @property
    def size(self):
        return tuple(reversed(self.operands[0].size))

    @property
    def rows(self):
        return self.operands[0].columns

    @property
    def columns(self):
        return self.operands[0].rows

    @property
    def operand(self):
        return self.operands[0]

    @operand.setter
    def operand(self, value):
        self.operands[0] = value

    @property
    def bandwidth(self):
        # Bandwidth is not switched/reversed here: The inverse of a lower
        # triangular matrix is still lower triangular.
        bands = self.operands[0].bandwidth
        return utils._inverse_bandwidth(bands, self.size)

    def __str__(self):
        return "{0}{1}".format(self.operands[0], self.name)

    def to_julia_expression(self, recommended=False, strip_inverse=False):
        if strip_inverse:
            return self.operands[0].to_julia_expression(recommended)
        else:
            return "inv({0})".format(self.operands[0].to_julia_expression(recommended))

    def to_cpp_expression(self, lib, recommended=False, strip_inverse=False):
        op = self.operands[0].to_cpp_expression(lib, recommended)
        if strip_inverse:
            return op
        return {
            CppLibrary.Blaze : "blaze::inv({0})",
            CppLibrary.Eigen : "({0}).inverse()",
            CppLibrary.Armadillo : "({0}).i()"
        }.get(lib).format(op)

class InverseTranspose(Operator):
    """docstring for InverseTranspose"""

    name = "^-T"
    arity = matchpy.Arity.unary
    infix = False

    def __init__(self,  *operands, variable_name=None):
        super(InverseTranspose, self).__init__(*operands, variable_name=variable_name)

    @property
    def size(self):
        return self.operands[0].size

    @property
    def rows(self):
        return self.operands[0].rows

    @property
    def columns(self):
        return self.operands[0].columns

    @property
    def operand(self):
        return self.operands[0]

    @operand.setter
    def operand(self, value):
        self.operands[0] = value

    @property
    def bandwidth(self):
        bands = tuple(reversed(self.operands[0].bandwidth))
        return utils._inverse_bandwidth(bands, self.size)

    def __str__(self):
        return "{0}{1}".format(self.operands[0], self.name)

    def to_julia_expression(self, recommended=False, strip_inverse=False):
        if strip_inverse:
            return "{0}'".format(self.operands[0].to_julia_expression(recommended))
        else:
            return "inv({0}')".format(self.operands[0].to_julia_expression(recommended))

    def to_cpp_expression(self, lib, recommended=False, strip_inverse=False):
        op = self.operands[0].to_cpp_expression(lib, recommended)
        if strip_inverse:
            return op
        return {
            CppLibrary.Blaze : "blaze::inv(blaze::trans({0}))",
            CppLibrary.Eigen : "({0}).transpose().inverse()",
            CppLibrary.Armadillo : "({0}).t().i()"
        }.get(lib).format(op)

class InverseConjugate(Operator):
    """docstring for InverseConjugate"""

    name = "^-C"
    arity = matchpy.Arity.unary
    infix = False

    def __init__(self,  *operands, variable_name=None):
        super(InverseConjugate, self).__init__(*operands, variable_name=variable_name)

    @property
    def size(self):
        return tuple(reversed(self.operands[0].size))

    @property
    def rows(self):
        return self.operands[0].columns

    @property
    def columns(self):
        return self.operands[0].rows

    @property
    def operand(self):
        return self.operands[0]

    @operand.setter
    def operand(self, value):
        self.operands[0] = value

    @property
    def bandwidth(self):
        bands = self.operands[0].bandwidth
        return utils._inverse_bandwidth(bands, self.size)

    def __str__(self):
         return "{0}{1}".format(self.operands[0], self.name)

class InverseConjugateTranspose(Operator):
    """docstring for InverseConjugateTranspose"""

    name = "^-H"
    arity = matchpy.Arity.unary
    infix = False

    def __init__(self,  *operands, variable_name=None):
        super(InverseConjugateTranspose, self).__init__(*operands, variable_name=variable_name)

    @property
    def size(self):
        return self.operands[0].size

    @property
    def rows(self):
        return self.operands[0].rows

    @property
    def columns(self):
        return self.operands[0].columns

    @property
    def operand(self):
        return self.operands[0]

    @operand.setter
    def operand(self, value):
        self.operands[0] = value

    @property
    def bandwidth(self):
        bands = tuple(reversed(self.operands[0].bandwidth))
        return utils._inverse_bandwidth(bands, self.size)

    def __str__(self):
         return "{0}{1}".format(self.operands[0], self.name)

######################
# symbols
#


class Scalar(Symbol):
    """docstring for Scalar"""
    def __init__(self, name, indices=set()):
        super(Scalar, self).__init__(name, (1, 1), indices)


class Vector(Symbol):
    """docstring for Vector"""
    def __init__(self, name, size, indices=set()):
        super(Vector, self).__init__(name, size, indices)
        if (size[0]==1 and size[1]==1) or (size[0]!=1 and size[1]!=1):
            raise ValueError("Symbol with size {} is not a vector.".format(size))

    def __repr__(self):
        if self.size[0] == 1:
            return 'Vector({!r}, cols={})'.format(self.name, self.size[1])
        else:
            return 'Vector({!r}, rows={})'.format(self.name, self.size[0])


class Matrix(Symbol):
    """docstring for Matrix"""
    def __init__(self, name, size, indices=set()):
        super(Matrix, self).__init__(name, size, indices)

    def __repr__(self):
        return 'Matrix({!r}, size={})'.format(self.name, self.size)


######################
# constant symbols
#


class ConstantMatrix(Matrix, Constant):
    """docstring for ConstantMatrix"""
    def __init__(self, name, size):
        super(ConstantMatrix, self).__init__(name, size)


class IdentityMatrix(ConstantMatrix):
    """docstring for IdentityMatrix"""
    def __init__(self, rows, columns):
        super(IdentityMatrix, self).__init__("I", (rows, columns))
        self.set_property(properties.IDENTITY)
        self.set_property(properties.DIAGONAL)

    def __repr__(self):
        return "{0}({1}, {2})".format(self.__class__.__name__, *self.size)

    def to_julia_expression(self):
        return "eye({0}, {1}, {2})".format(config.data_type_string, *self.size)


class ZeroMatrix(ConstantMatrix):
    """docstring for ZeroMatrix"""
    def __init__(self, rows, columns):
        super(ZeroMatrix, self).__init__("0", (rows, columns))
        self.set_property(properties.ZERO)

    def __repr__(self):
        return "{0}({1}, {2})".format(self.__class__.__name__, *self.size)

    def to_julia_expression(self):
        return "zeros({0}, {1}, {2})".format(config.data_type_string, *self.size)

class ConstantScalar(Scalar, Constant):
    """docstring for ConstantScalar"""
    def __init__(self, value):
        super(ConstantScalar, self).__init__(str(value))
        self.value = value

    def __repr__(self):
        return "{0}({1})".format(self.__class__.__name__, self.value)

    def to_julia_expression(self):
        return "{0}".format(self.value)


class Index(object):
    def __init__(self, name, value):
        self.name = name
        self.value = value
        self._hash = hash((name, value))

    @staticmethod
    def range(indices):
        return functools.reduce(operator.mul, [index.value for index in indices])

    def __eq__(self, other):
        return self.value == other.value and self.name == other.name

    def __hash__(self):
        return self._hash

    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name


class BlockedExpression(Expression):
    def __init__(self, nd_array):
        # TODO check consistency here? probably yes
        Expression.__init__(self, *nd_array)
        self.check_validity()

    @property
    def size(self):
        rows = sum(row[0].size[0] for row in self.operands)
        cols = sum(block.size[1] for block in self.operands[0])
        return (rows, cols)

    @property
    def rows(self):
        return sum(row[0].rows for row in self.operands)

    @property
    def columns(self):
        return sum(block.columns for block in self.operands[0])

    def check_validity(self):
        # Check if rows have the same number of cells.
        row_lengths = [len(row) for row in self.operands]
        if not row_lengths.count(row_lengths[0]) == len(row_lengths):
            raise ExpressionError("Rows in " + repr(self) + " do not have the same number of cells.")

        # Check if rows have consistent sizes.
        for i, row in enumerate(self.operands):
            row_sizes = [block.size[0] for block in row]
            if not row_sizes.count(row_sizes[0]) == len(row_sizes):
                raise ExpressionError("Sizes in row " + str(i) + " (" + str(row) + ") in " + repr(self) + " do not match.")

        # Check if columns have consistent sizes.
        operands_transposed = list(zip(*self.operands))
        for i, column in enumerate(operands_transposed):
            column_sizes = [block.size[1] for block in column]
            if not column_sizes.count(column_sizes[0]) == len(column_sizes):
                raise ExpressionError("Sizes in colum " + str(i) + " (" + str(list(column)) + ") in " + repr(self) + " do not match.")

    def set_child(self, i, expr):
        #pointer = self.operands
        #for pos_i in range(len(i)-1):
            #pointer = pointer[ i[pos_i] ]
        #pointer[i[-1]] = expr
        row, col = i // len(self.operands[0]), i % len(self.operands[0])
        self.operands[row][col] = expr

    # only works for flattening matrices
    def flatten_operands(self):
        return list(itertools.chain.from_iterable(self.operands))

    def _cleanup(self):
        # TODO: should reassign back to operands
        for ch in self.flatten_operands():
            ch._cleanup()
        return self

    def __iter__(self):
        yield from self.operands

    def match(self, ctx):
        # nothing to match, failure
        if len(ctx.stack_expr) == 0:
            return None
        # pop the expression to match
        expr = ctx.stack_expr.pop()
        if self.__class__ == expr.__class__ and \
                self.shape == expr.shape:
            # This allows the use of WildcardPlus and WildcardStar for rows or full blocked expressions
            # In principle, no BlockedExpression would appear in a pattern, would it?
            patt_seq = Sequence(*self.flatten_operands())
            expr_seq = Sequence(*expr.flatten_operands())
            ctx.stack_expr.append(expr_seq)
            for m in patt_seq.match(ctx):
                yield m

    def iterate_preorder(self):
        yield self
        for child in self.flatten_operands():
            yield from child.iterate_preorder()

    # only for matrices (2D blocked expressions)
    def _preorder_position(self, parent=(None, None)):
        yield (id(self), parent)
        for i in range(len(self.operands)):
            for j in range(len(self.operands[0])):
                yield from self.operands[i][j]._preorder_position((self, (i,j)))

    def _postorder_stack(self, parent=(None,None)):
        return [ self, parent, [ch._postorder_stack((self, i)) for i, ch in enumerate(self.flatten_operands())] ]

    def __eq__(self, other):
        return self.__class__ == other.__class__ and \
                self.operands == self.operands

    def __repr__(self):
        #return "[ %s ]" % ("; ".join([ ", ".join([ cell for cell in row ]) for row in self.operands ]))
        return repr(self.operands)

    def __str__(self):
        return str(self.operands)

if __name__ == "__main__":

    n = 10

    A = Matrix("A", (10, 10))
    A.set_property(properties.LOWER_TRIANGULAR)
    print(A.has_property(properties.LOWER_TRIANGULAR))
    B = Matrix("B", (10, 10))
    C = Matrix("C", (10, 10))
    CC = ConstantMatrix("CC", (10, 10))

    I = IdentityMatrix(10, 10)

    fI = matchpy.freeze(I)
    print(repr(fI))
    print(type(fI))

    two = ConstantScalar(2)

    ftwo = matchpy.freeze(two)
    print(repr(ftwo))
    print(type(ftwo))
    print(Plus)

    # ftwo.value = 4

    # expr = Times(A, B)
    expr = matchpy.freeze(Plus(B, Times(B, C)))
    # expr = two
    print(str(expr))
    print(repr(expr))
    print(repr(type(expr)))
    print(repr(type(expr.operands[0])))

    expr = Times(A, Times(B, Inverse(C)))
    expr = matchpy.freeze(Times(A, Times(B, Inverse(C))))

    print(expr.has_property(properties.LOWER_TRIANGULAR))
    print(expr.bandwidth)
    # print(Times(A, B))
    # print(Times(A))
    # print(Identity(A))
    # print(str(I))
    # print(repr(I))