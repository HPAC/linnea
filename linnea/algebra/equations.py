from .. import temporaries
from ..derivation.partitioning import _propagate_partitioning, apply_partitioning

from . import expression as ae
from . import transformations as at
from . import representations as ar

from .properties import Property, binary_implications
from .validity import check_validity
from .consistency import check_consistency

import copy
import matchpy

class UnknownSymbolType(Exception):
    pass

class InvalidInput(Exception):
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

    def set_equivalent(self, equations_before):
        """Applies temporaries.set_equivalent() to all equations.
 
        Args:
            equations_before (Equations)       
        """
        for n, equation in enumerate(self.equations):
            if equation.rhs != equations_before[n].rhs:
                temporaries.set_equivalent(equations_before[n].rhs, equation.rhs)

    def remove_identities(self):
        """Removes equations where both sides are temporaries.

        There are three different cases.

        If both sides of an equation are the same, it can simply be removed.

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

        Furthermore, it is assumed that there are no cases such as
        tmp2 = tmp3
        tmp1 = tmp2

        Returns:
            Equations: self with equations removed.
        """
        remove = []
        replace_eqns = []
        equations = list(self.equations)
        for n, equation in enumerate(equations):
            if equation.rhs.name in temporaries._equivalent_expressions:
                if equation.lhs.name in temporaries._equivalent_expressions:
                    remove.append(n)
                    if equation.lhs != equation.rhs:
                        replace_eqns.append(equation)
                else:
                    # TODO it would be better to only do this if the operand is known to be an intermediate.
                    replace_eqns.append(equation)

        for idx in reversed(remove):
            del equations[idx]

        for replace_eqn in replace_eqns:
            rule = matchpy.ReplacementRule(matchpy.Pattern(replace_eqn.lhs), lambda: replace_eqn.rhs)
            equations_replaced = []
            for equation in equations:
                equations_replaced.append(ae.Equal(equation.lhs, matchpy.replace_all(equation.rhs, (rule,))))
            equations = equations_replaced

        # TODO do we want to manipulate the table of temporaries here?
        return Equations(*equations)

    def replace_intermediates(self):
        candidates = []
        for equation in self.equations:
            if (equation.rhs.name in temporaries._equivalent_expressions
                and not equation.lhs.name in temporaries._equivalent_expressions):
                candidates.append(equation)

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
        input, _ = self.input_output()
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


    def check_validity(self):
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
        
        return all(check_validity(equation) for equation in self.equations)


    def check_consistency(self):
        return all(check_consistency(equation) for equation in self.equations)        


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
