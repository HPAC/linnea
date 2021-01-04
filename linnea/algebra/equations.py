from .. import temporaries
from ..derivation.partitioning import _propagate_partitioning, apply_partitioning

from . import expression as ae
from . import transformations as at
from . import representations as ar

from .properties import Property, binary_implications
from .validity import check_validity

import copy
import matchpy

class UnknownSymbolType(Exception):
    pass

class InvalidInput(Exception):
    pass


class Replacer(object):
    """An alternative for lambda functions in ReplacementRule objects.

    This object can be used in ReplacementRule objects instead of constant
    lambda functions. It avoids problems that can occur if the contents of
    variables used in a lambda function changes in the current scope.
    """
    def __init__(self, replacement):
        self.replacement = replacement

    def __call__(self, **k):
        return self.replacement


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
            equation = ar.to_normalform(equation)
            equations.append(equation)
        return Equations(*equations)

    def simplify(self):
        return Equations(*map(at.simplify, self.equations))

    def to_SOP(self):
        return Equations(*map(ar.to_SOP, self.equations))

    def to_normalform(self):
        return Equations(*map(ar.to_normalform, self.equations))

    def get_data(self):
        """Counts the number of floating point values for calculating intensity.

        
        """

        operands_lhs = set()
        operands_rhs = set()
        for equation in self.equations:
            operands_lhs.update(op for op, _ in equation.lhs.preorder_iter() if isinstance(op, ae.Symbol))
            operands_rhs.update(op for op, _ in equation.rhs.preorder_iter() if isinstance(op, ae.Symbol))
        
        all_data = 0
        for operands in (operands_lhs, operands_rhs):
            data = 0
            for operand in operands:
                rows, cols = operand.size
                if not operand.has_property(Property.CONSTANT):
                    if operand.has_property(Property.DIAGONAL):
                        data += min(rows, cols)
                    elif operand.has_property(Property.SYMMETRIC) or operand.has_property(Property.TRIANGULAR):
                        data += rows * cols/2
                    else:
                        data += rows * cols
            all_data += data

        return all_data

    def remove_identities(self):
        """Removes equations of the form tmpX = tmpY.

        There are three different cases:
        
        If both sides of the equation are the same temporary, the equation can
        simply be removed.

        If both sides are different temporaries, all occurrences of the
        temporary on the left-hand side in the other equations are replaced by
        the temporary on the right-hand side. This is necessary to ensure that
        code generation still works. This case can happen when setting
        equivalent expressions does not work perfectly.

        If the operand on the left-hand side is an intermediate, and the one on
        the right hand side is a temporary, all occurrences of the intermediate
        will be replace with the temporary.

        Replacing the temporary in the other equations only works if the
        temporary that is replaced is not part of any computation yet.

        Replacing temporaries even works if there are sequences of assignments
        such as:
        tmp2 = tmp3
        tmp1 = tmp2

        Returns:
            Equations: self without identities, and with replaced temporaries.
        """

        new_equations = []
        mapping = dict()
        for equation in self.equations:
            if temporaries.is_temporary(equation.rhs):
                if equation.lhs != equation.rhs:
                    try:
                        op = mapping[equation.rhs]
                    except KeyError:
                        mapping[equation.lhs] = equation.rhs
                    else:
                        mapping[equation.lhs] = op
                    if temporaries.is_temporary(equation.lhs):
                        continue
                else:
                    continue

            if mapping:
                rules = [matchpy.ReplacementRule(matchpy.Pattern(rhs), Replacer(lhs)) for rhs, lhs in mapping.items()]
                new_equations.append(ae.Equal(equation.lhs, matchpy.replace_all(equation.rhs, rules)))
            else:
                new_equations.append(equation)

        return Equations(*new_equations)

    def remove_explicit_transposition(self, eqn_idx=None):
        """Removes equations of the form tmpX = tmpY^T.

        With common subexpression elimination and the application of tranposed
        operands, it is possible to reach assignments of the form tmpX = tmpY^T.
        To avoid that an explicit transposition is computed, this function
        replaces tmpX in all subsequent assignments with tmpY^T.

        Args:
            eqn_idx (int, optional): Only test if self.equations[eqn_idx] is an
                explicit transposition.

        Returns:
            Equations: self with explicit transposition removed.
        """              

        indices = None
        if eqn_idx:
            indices = [eqn_idx]
        else:
            indices = range(len(self.equations))

        replacement_rules = []
        remove = set()
        for idx in indices:
            equation = self.equations[idx]
            if temporaries.is_temporary(equation.lhs) and isinstance(equation.rhs, ae.Transpose) and isinstance(equation.rhs.operand, ae.Symbol):
                if equation.lhs != equation.rhs.operand:
                    # In they are equal, the operand is symmetric, which means
                    # that this assignment is an identity and can be removed.
                    replacement_rules.append(matchpy.ReplacementRule(matchpy.Pattern(equation.lhs), Replacer(equation.rhs)))
                remove.add(idx)

        if replacement_rules:
            new_equations = []
            for i, equation in enumerate(self.equations):
                if not i in remove:
                    if i < min(remove):
                        # It is sufficient to start replacing at the smallest index in remove.
                        new_equations.append(equation)
                    else:
                        new_equations.append(ae.Equal(equation.lhs, matchpy.replace_all(equation.rhs, replacement_rules)))
            return Equations(*new_equations)
        else:
            return self

    def process_next(self):
        for eqn_idx, equation in enumerate(self.equations):
            if not isinstance(equation.rhs, ae.Symbol):
                return equation, eqn_idx
        return None, -1

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
                if isinstance(node, ae.Symbol) and not node.has_property(Property.CONSTANT):
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
            if eqn.lhs.has_property(Property.AUXILIARY):
                rule = (matchpy.Pattern(eqn.lhs), lambda **_: eqn.rhs)
                replacement_rules.append(rule)
                remove.append(eqn)

        for eqn in remove:
            self.equations.remove(eqn)

        self.equations = [matchpy.replace_all(equation, replacement_rules) for equation in self.equations]


    def infer_lhs_properties(self):
        for equation in self.equations:
            operand = equation.lhs
            if isinstance(operand, ae.Symbol):
                # Checking if the operand is a symbol is necessary because when
                # this function is called, correctnes has not been checked yet.
                for property in Property:
                    if equation.rhs.has_property(property):
                        operand.set_property(property)


    def infer_missing_properties(self):
        """Infers missing properties of operands.

        Infers properties of input operands that follow from combinations of
        other properties. Those properties might be missing if the user did not
        use the most specific properties.

        In addition, this function adds the properties SQUARE, ROW_PANEL, and
        COLUMN_PANEL. The more elegant solution would be to set those properties
        in Matrix.__init__(). However, this is currently not possible because
        the properties of matrices are used in the kernel description to
        generate the constraints for variables. There are many kernels that do
        not require a specific shape, but since we are using actual Matrix
        objects for the kernel description, they do have exactly one of those
        three shapes, leading to unwanted constraints.
        To fix this, the way KernelDescription objects work needs to be changed.
        """
        input, output = self.input_output()
        for operand in input:
            
            rows, columns = operand.size
            if rows < columns:
                operand.set_property(Property.ROW_PANEL)
            elif rows > columns:
                operand.set_property(Property.COLUMN_PANEL)
            elif rows != 1: # Scalars must not be square.
                operand.set_property(Property.SQUARE)

            new_property = True
            while new_property:
                for property_set, property in binary_implications.items():
                    if property not in operand.properties and property_set <= operand.properties:
                        operand.properties.add(property)
                        new_property = True
                        break
                else:
                    new_property = False

        for operand in output:
            if isinstance(operand, ae.Symbol):
                # Checking if the operand is a symbol is necessary because when
                # this function is called, correctnes has not been checked yet.
                rows, columns = operand.size
                if rows < columns:
                    operand.set_property(Property.ROW_PANEL)
                elif rows > columns:
                    operand.set_property(Property.COLUMN_PANEL)
                elif rows != 1: # Scalars must not be square.
                    operand.set_property(Property.SQUARE)

    def check_validity(self, dependent_dimensions=False):
        """Checks if the object is valid.

        In addition to checking that the individual equations are well-formed,
        this function also checks that the sequence of equations is valid.
        All operands are treated as unique mathematical objects, that this,
        they only represent one value. As a result, 1) operands cannot be
        redefined, and 2) input operands cannot be overwritten. Since at
        present, Linnea does not use any Property or labels to annotate what
        is input and what is output, the two requirements are checked by
        analyzing the dataflow: Requirement 2) is checked by testing that the
        first use of an operands takes place after its definition (if there is
        a definition; otherwise, the operand is an input operand).

        Examples of invalid sequences of equations:
        X = X A

        Y = X A
        X = B C

        X = A
        X = B

        Args:
            dependent_dimensions (bool, optinal): If True, disables some checks
                that are not applicable if the input is only used to compute
                dependent dimensions.

        Returns:
            bool: True if the equations are valid, raises an exception
                otherwise.

        Raises:
            InvalidInput: If any of the two requirements above is violated.
            SizeMismatch: If operand sizes in any expression do not match.

        """
        first_use = dict()
        definition = dict()
        for n, equation in enumerate(self.equations):
            # Test that there is only on definition.
            if not equation.lhs in definition:
                definition[equation.lhs] = n
            else:
                raise InvalidInput("Operands can only be defined once: {}".format(equation.lhs.name)) 

            for expr, _ in equation.rhs.preorder_iter():
                if isinstance(expr, ae.Symbol) and not isinstance(expr, ae.Constant):
                    # add expr: n to first_use if expr is not in dict
                    if not expr in first_use:
                        first_use[expr] = n
        
        # Test that first use happens after definition.         
        for operand, definition_line in definition.items():
            try:
                first_use_line = first_use[operand]
            except KeyError:
                pass
            else:
                if first_use_line <= definition_line:
                    raise InvalidInput("Input operands can not be overwritten: {}".format(operand.name))
        return all(check_validity(equation, dependent_dimensions) for equation in self.equations)

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

    def to_matlab_expression(self):
        return "\n".join([equation.to_matlab_expression() for equation in self.equations])

    def to_julia_expression(self):
        return "\n".join([equation.to_julia_expression() for equation in self.equations])

    def to_cpp_expression(self, lib):
        return "\n".join([equation.to_cpp_expression(lib) for equation in self.equations])

    def to_tex(self):
        strings = []
        for equation in self.equations:
            strings.append(" ".join([equation.lhs.to_tex(), "&:=", equation.rhs.to_tex()]))
        return "\\begin{{align}}\n{}\n\\end{{align}}".format("\\\\\n".join(strings))

    def input_output(self):
        input = []
        output = []
        seen_before = set()
        for equation in self.equations:
            output.append(equation.lhs)
            for expr, _ in equation.rhs.preorder_iter():
                if isinstance(expr, ae.Symbol) and not isinstance(expr, ae.Constant) and not expr in seen_before and not expr in output:
                    input.append(expr)
                    seen_before.add(expr)
        return input, output
