from . import expression as ae
from .properties import Property as properties

from . import transformations as at
from . import representations as ar

from ..derivation.partitioning import _propagate_partitioning, apply_partitioning

from .. import temporaries

import copy

import matchpy

class UnknownSymbolType(Exception):
    pass

class Equations():

    _counter = 0

    def __init__(self, *equations):
        self.equations = tuple(equations)

    def __iter__(self):
        return self.equations.__iter__()

    def __repr__(self):
        return "\n".join([repr(equation) for equation in self.equations])
        
    def __str__(self):
        return "\n".join([str(equation) for equation in self.equations])

    def __eq__(self, other):
        return hash(self) == hash(other)

    def __lt__(self, other):
        return self.equations < other.equations

    def __getitem__(self, key):
        return self.equations[key]

    def __len__(self):
        return len(self.equations)

    def __hash__(self):
        return hash(tuple(self.equations))

    def set(self, position, value):
        equations = list(self.equations)
        equations[position] = value
        return Equations(*equations)

    def insert(self, position, value):
        equations = list(self.equations)
        equations.insert(position, value)
        return Equations(*equations)

    def replace_all(self, rules):
        equations = []
        for equation in self.equations:
            equation = matchpy.replace_all(equation, rules)
            equation = ar.to_SOP(at.simplify(ar.to_SOP(equation)))
            equations.append(equation)
        return Equations(*equations)

    def simplify(self):
        return Equations(*[at.simplify(equation) for equation in self.equations])

    def to_SOP(self):
        return Equations(*[ar.to_SOP(equation) for equation in self.equations])

    def to_normalform(self):
        return Equations(*[ar.to_SOP(at.simplify(ar.to_SOP(equation))) for equation in self.equations])

    def set_equivalent(self, equations_before):
        """Applies temporaries.set_equivalent() to all equations.
 
        Args:
            equations_before (Equations)       
        """
        for n, equation in enumerate(self.equations):
            if equation.rhs != equations_before[n].rhs:
                temporaries.set_equivalent(equations_before[n].rhs, equation.rhs)

    def metric(self):
        # TODO how to compute the metric of multiple equations?
        sum = [0, 0]
        for equation in self.equations:
            m = equation.rhs.metric()
            sum[0] += m[0]
            sum[1] += m[1]
        return sum

    def remove_identities(self):
        """Removes equations where both sides are temporaries.

        There are two different cases. If both sides of an equation are the
        same, it can simply be removed.

        If both sides are different, all occurrences of the temporary on the
        left-hand side in the other equations are replaced by the temporary on
        the right-hand side. This is necessary to ensure that code generation
        still works. This case can happen when setting equivalent expressions
        does not work perfectly.

        Replacing the temporary in the other equations only works if the
        temporary that is replaced is not part of any computation yet.

        Furthermore, it is assumed that there are no cases such as
        tmp2 = tmp3
        tmp1 = tmp2

        Returns:
            Equations: self with equations removed.
        """
        remove = []
        replace_eqns = []
        equations = list(self.equations)
        for n, equation in enumerate(self.equations):
            if (equation.lhs.name in temporaries._equivalent_expressions 
                    and equation.rhs.name in temporaries._equivalent_expressions):
                remove.append(n)
                if equation.lhs != equation.rhs:
                    replace_eqns.append(equation)

        for idx in reversed(remove):
            del equations[idx]

        for replace_eqn in replace_eqns:
            rule = matchpy.ReplacementRule(matchpy.Pattern(replace_eqn.lhs), lambda: replace_eqn.rhs)
            equations = [matchpy.replace_all(equation, (rule,)) for equation in equations]

        # TODO do we want to manipulate the table of temporaries here?
        return Equations(*equations)

    def apply_partitioning(self):
        change = True
        # The partitioning can not be propagated per equation. If the same
        # operand shows up in multiple equations, that creates a dependency
        # in both directions between those equations. Thus, the partitioning
        # has to be propagated in all equations simultanously.
        while change:
            change = False
            for equation in self.equations:
                change = change or _propagate_partitioning(equation)
                
        self.equations = [apply_partitioning(equation) for equation in self.equations]
        equations = []
        for equation in self.equations:
            if isinstance(equation, ae.BlockedExpression):
                equations.extend(equation.flatten_children())
            else:
                equations.append(equation)
        self.equations = equations

    def resolve_dependencies(self):
        """Resolves dependencies between the equations.

        This function ensures that the same symbols does not represent different
        expressions. This is the case when variables are overwritten.
        Example: X = X A is transformed into X0 = X A
        """
        seen_before = set()
        replacements = dict()

        equations = []
        for equation in self.equations:
            # look at rhs
            for node, pos in equation.rhs.preorder_iter():
                if isinstance(node, ae.Symbol) and not node.has_property(properties.CONSTANT):
                    seen_before.add(node.name)
                    if node.name in replacements:
                        equation.rhs.set_successor(pos, replacements[node.name])
            # look at lhs
            lhs = equation.lhs
            if lhs.name in seen_before:
                new_symbol = self._copy_symbol(lhs)
                replacements[lhs.name] = new_symbol
                equation.lhs = new_symbol
            else:
                seen_before.add(lhs.name)

    def replace_auxiliaries(self):
        """Replaces auxiliaries.

        If there is an equation X = rhs where X has the property "auxiliary", the
        equation is removed and all occurences of X in other equations are
        replaced with rhs.
        """
        replacement_rules = []
        remove = []
        # IMPORTANT: eqn can NOT be renamed to equation. Otherwise,
        # "equation.rhs in the lambda function refers to the wrong equation.
        for eqn in self.equations:
            if eqn.lhs.has_property(properties.AUXILIARY):
                rule = (matchpy.Pattern(eqn.lhs), lambda **_: eqn.rhs)
                replacement_rules.append(rule)
                remove.append(eqn)

        for eqn in remove:
            self.equations.remove(eqn)

        self.equations = [matchpy.replace_all(equation, replacement_rules) for equation in self.equations]


    def infer_lhs_properties(self):
        output_operands = set()
        for equation in self.equations:
            operand = equation.lhs
            for property in properties:
                if equation.rhs.has_property(property):
                    operand.set_property(property)


    def _copy_symbol(self, symbol):
        """Creates a copy of symbol with a unique name"""
        new_symbol = None
        name = "".join([symbol.name, str(Equations._counter)])
        Equations._counter += 1
        if isinstance(symbol, ae.Scalar):
            new_symbol = ae.Scalar(name, symbol.size)
        if isinstance(symbol, ae.Vector):
            new_symbol = ae.Vector(name, symbol.size)
        elif isinstance(symbol, ae.Matrix):
            new_symbol = ae.Matrix(name, symbol.size)
        else:
            print(symbol)
            raise UnknownSymbolType()

        return new_symbol

    def to_dot_files(self):
        for equation in self.equations:
            equation.to_dot_file()

    def to_matlab_expression(self, recommended=False):
        return "\n".join([equation.to_matlab_expression(recommended) for equation in self.equations])

    def to_julia_expression(self, recommended=False):
        return "\n".join([equation.to_julia_expression(recommended) for equation in self.equations])

    def to_cpp_expression(self, lib, recommended=False):
        return "\n".join([equation.to_cpp_expression(lib, recommended) for equation in self.equations])

    def input_output(self):
        input = []
        output = []
        seen_before = set()
        for equation in self.equations:
            output.append(equation.lhs)
            for expr, _ in equation.rhs.preorder_iter():
                if isinstance(expr, ae.Symbol) and not isinstance(expr, ae.Constant) and not expr in seen_before:
                    input.append(expr)
                    seen_before.add(expr)
        return input, output
