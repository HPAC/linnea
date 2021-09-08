from ..algebra.expression import Symbol, Times, \
                                 Inverse, InverseTranspose, Transpose, \
                                 LinSolveL, LinSolveR, \
                                 Matrix, Vector, Scalar
from ..algebra.equations import Equations
from ..algebra.properties import Property
from ..utils import PropertyConstraint, is_inverse
from .. import kernels

from .. import config

from .memory import memory as memory_module

import copy
import math
import textwrap
import os
import matchpy
import pkg_resources
import itertools

class Algorithm():
    """Represents an Algorithm and translates it to code.

    This class represents one algorithm and offers functions to translate it to
    code.

    Attributes:
        name (string): The name of the algorithm.
        initial_equations (Equations): The input equations.
        final_equations (Equations): The equations left at the end of the
            dervation. They only consist of quations of the form
            Symbol = Symbol. Is used for the mapping of temporaries to output
            operands.
        kernels_and_equations (list): List of MatchedKernel and Equation
            objects, representing the algorithm (and its generation steps).
        cost (int): The cost of the entire algorithm.
        data (int): The amount of data (number of scalars) of the input
            equations.
        intensity (float): The computational intensity of this algorithm.
        liveness (dict): Stores liveness information. For each operand, it
            contains a tuple of to integers. The first integer is the line
            number (corresponds to a MatchedKernel object in matched_kernels)
            where the operand is define. If it is None, the operand is an input
            operand. The second integer is the line number where the operand is
            used for the last time. If it is None, the operand is an ouput
            operand. Is written when accessing the code property.
        docstring (string): Docstring for the generated code. Is written when
            accessing the code property.
    """
    def __init__(self, name, initial_equations, final_equations, kernels_and_equations, cost):

        # super(Algorithm, self).__init__()
        self.name = name
        self.initial_equations = initial_equations
        self.final_equations = final_equations
        self.matched_kernels = [keq for keq in kernels_and_equations if isinstance(keq, MatchedKernel)]
        self.kernels_and_equations = kernels_and_equations
        self.cost = cost
        self.data = self.initial_equations.get_data()
        self.intensity = self.cost/self.data

        self.liveness = None

        self._experiment_input = None
        self._experiment_output = None

        self._code = None
        self.docstring = None

    def __eq__(self, other):
        return self.matched_kernels == other.matched_kernels

    def __hash__(self):
        return hash(tuple(self.matched_kernels))

    def _symbols(self, expr):
        """Returns all the symbols in expr.
    
        """
        # TODO this should go into some utils module
        if isinstance(expr, Symbol):
            yield expr
        else:
            for operand in expr.operands:    
                yield from self._symbols(operand)

    @property
    def code(self):
        """str: The generated code."""
        if self._code is None:
            self._code = self._generate_code()
        return self._code

    def code_as_function(self, experiment=False):
        """Returns the generated code as a function.
        
        Args:
            experiment (bool, optional): If True, the generated function
                contains a timer. Defaults to False.

        Returns:
            string: The generated code.
        """
        if config.language == config.Language.Julia:
            code = self.code # this cannot be done later because it sets self.docstring, self._experiment_input, and self._experiment_output
            if experiment:
                template = get_template("algorithm_experiments.jl", config.language)
                algorithm_str = template.format(self.name, self._experiment_input, textwrap.indent(code, "    "), self._experiment_output)
            else:
                template = get_template("algorithm.jl", config.language)
                algorithm_str = template.format(self.docstring, self.name, self._experiment_input, textwrap.indent(code, "    "), self._experiment_output)
        else:
            raise config.LanguageOptionNotImplemented()
        return algorithm_str

    def liveness_analysis(self):
        """Performs a liveness analysis.

        Constructs and sets the self.liveness dictionary by traversing the
        matched kernels in reversed order. If a symbols shows up on the RHS of
        an equation for the first time, that is its last use. If a symbols shows
        up on the LHS of an equation, that is its definition (this requires that
        the operations of matched_kernels are in static single assignment form,
        which is the case because they are purely symbolic, and a symbol
        represents a unique operand).
        """
        liveness = dict()
        live_symbols = set()
        for line_number, matched_kernel in reversed(list(enumerate(self.matched_kernels))):
            rhs_symbols = set(self._symbols(matched_kernel.operation.rhs))
            for rhs_symbol in rhs_symbols:
                if rhs_symbol not in live_symbols:
                    liveness.setdefault(rhs_symbol.name, [None, None])[1] = line_number
            
            live_symbols.update(rhs_symbols)

            lhs_symbols = set(self._symbols(matched_kernel.operation.lhs))
            for lhs_symbol in lhs_symbols:
                liveness.setdefault(lhs_symbol.name, [None, None])[0] = line_number

        for equation in self.final_equations:
            try:
                liveness[equation.rhs.name][1] = None
            except KeyError:
                # If equation.rhs.name is not in liveness at this point, it must
                # be an input operand that is not used in any kernel call. As a
                # result, this operand is both input and output operand.
                liveness[equation.rhs.name] = [None, None]

        return liveness

    def _generate_code(self):
        """Translated the algorithm to code.

        The language depends on what is specified in config.
        
        Returns:
            string: Code.
        """


        """
        TODO this is inconsistent. memory and known_lines are attributes,
        matched_kernel and line_number is passed to memory.add_operation
        TODO it doesn't make a lot of sense that memory is part of the object,
        same for known_lines
        """

        memory = memory_module.Memory(self.initial_equations)
        self.liveness = self.liveness_analysis()
        
        known_lines = set()
        code_list = []

        input_operands, output_operands = self.initial_equations.input_output()
        self._experiment_input = ", ".join(["{}::{}".format(memory.lookup[operand.name].name, operand_type(operand)) for operand in input_operands])

        self.docstring = self._generate_docstring(memory)

        code_list.append("{0}cost: {1:.3g} FLOPs\n".format(config.comment, self.cost))

        if config.c:
            code_list.append("int info = 0;\n\n")

        for line_number, matched_kernel in enumerate(self.matched_kernels):
            code_list.append(self._matched_kernel_to_code(matched_kernel, line_number, memory, known_lines))

        code_list.append(config.comment)
        code_list.append(memory.content_string_with_format())
        code_list.append("\n")

        output_operand_mapping = {eqn.lhs: eqn.rhs for eqn in self.final_equations}
        conversion_to_full = memory.convert_to_full(output_operand_mapping.values())
        if conversion_to_full:
            code_list.append(conversion_to_full)
            code_list.append("\n")

            code_list.append(config.comment)
            code_list.append(memory.content_string_with_format())
            code_list.append("\n")

        code_list.append(textwrap.indent(str(self.final_equations), config.comment))
        self._experiment_output = ", ".join([memory.lookup[output_operand_mapping[operand].name].name for operand in output_operands])

        """Resetting the MemoryLocation ID counter is purely cosmetic. Since
        they only need to be unique per algorithm and are not used for anything
        else, it is fine to do this here."""
        memory_module.MemoryLocation.reset_ID_counter()
        return "".join(code_list)

    def _matched_kernel_to_code(self, matched_kernel, line_number, memory, known_lines):
        """Generates code for a single MatchedKernel objects.

        Args:
            matched_kernel (MatchedKernel)
            line_number (int): The line number of the operation represeted by
                matched_kernel. Is necessary to use the liveness dictionary.
            memory (Memory): Represents the current state of the memory.
            known_lines (set): Contains lines of code that may be generated
                multiple times even though they are only needed once, such as
                definitions of constant arguments.

        Returns:
            string: Code.

        Important:
            MatchedKernel objects can not be modified because they are used
            multiple times for the generation of algorithms.
            The memory object is modified by this function.

        """
        lines_list = []

        # print("#########")
        # print(matched_kernel.signature)
        # print(matched_kernel.operand_dict)
        # print(matched_kernel.kernel_io)

        mem_content = "".join([config.comment, memory.content_string_with_format(), "\n"])
        # mem_content = "".join([config.comment, memory.content_string(), "\n"])
        lines_list.append(mem_content)

        # TODO arguments could be a stack or something, removing processed arguments. Then we don't need late_arguments.
        arguments = copy.copy(matched_kernel.arguments)
        argument_mapping = dict()
        argument_pre_code = []

        signature = matched_kernel.signature.safe_substitute_copy(matched_kernel.other_replacements)

        if matched_kernel.pre_code:
            kernel_pre_code = matched_kernel.pre_code.safe_substitute_copy(matched_kernel.other_replacements)
        if matched_kernel.post_code:
            kernel_post_code = matched_kernel.post_code.safe_substitute_copy(matched_kernel.other_replacements)

        late_arguments = []
        for argument in arguments:
            if isinstance(argument, kernels.utils.general.StorageFormatArgument):
                late_arguments.append(argument)
                continue
            try:
                pre_code, arg_replacement = argument.get_replacement(matched_kernel.operand_dict, memory)
            except memory_module.OperandNotInMemory:
                """ This Argument object has to be used after
                memory.add_operation(). Happens only for StrideArgument objects
                if nothing is overwritten.
                TODO: If nothing is overwritten, a new memory location is
                allocated. Thus, the stride can be chosen. So does  this  make 
                sense?
                """
                late_arguments.append(argument)
            else:
                if pre_code:
                    argument_pre_code.append(pre_code)
                argument_mapping[argument.name] = arg_replacement   

        # print(matched_kernel.operation)
        mem_ops_before, mem_ops_after, operand_mapping, actual_storage_formats = memory.add_operation(matched_kernel.kernel_io, self.liveness, line_number)

        if mem_ops_before:
            mem_code_before = "".join([mem_op.code() for mem_op in mem_ops_before])
            lines_list.append(mem_code_before)

        # print(late_arguments)
        for argument in late_arguments:
            pre_code, arg_replacement = argument.get_replacement(matched_kernel.operand_dict, memory, actual_storage_formats)
            if pre_code:
                argument_pre_code.append(pre_code)
            argument_mapping[argument.name] = arg_replacement 

        # print(argument_mapping)
        # print(signature)

        signature.safe_substitute(operand_mapping)

        if matched_kernel.pre_code:
            kernel_pre_code.safe_substitute(operand_mapping)
        if matched_kernel.post_code:
            kernel_post_code.safe_substitute_copy(operand_mapping)

        if argument_pre_code:
            lines_list.extend(self.remove_duplicate_lines(argument_pre_code, known_lines))

        if matched_kernel.pre_code:
            lines_list.append(kernel_pre_code.safe_substitute_str(argument_mapping))

        lines_list.append(config.comment)
        lines_list.append(str(matched_kernel.operation))
        lines_list.append("\n")

        lines_list.append(signature.safe_substitute_str(argument_mapping))
        lines_list.append("\n")

        if matched_kernel.post_code:
            lines_list.append(kernel_post_code.safe_substitute_str(argument_mapping))

        if mem_ops_after:
            mem_code_after = "".join([mem_op.code() for mem_op in mem_ops_after])
            lines_list.append(mem_code_after)

        lines_list.append("\n")
        return "".join(lines_list)


    def _matched_kernel_to_generation_steps(self, matched_kernel):
        """Generates pseudocode for a signle MatchedKernel objects.

        Args:
            matched_kernel (MatchedKernel)

        Important:
            MatchedKernel objects can not be modified because they are used
            multiple times for the generation of algorithms.

        Returns:
            string: Pseudocode.
        """
        # TODO Perhaps this should be a method of MatchedKernel.
        return "{0:<30}# {1:.3g}".format(str(matched_kernel.operation), matched_kernel.cost)

    def generation_steps(self):
        """Generates a description of how the algorithm was found.
        
        Returns:
            string: Description of the generation steps.
        """
        code_list = []

        code_list.append("# cost {:.3g}".format(self.cost))

        for kernel_or_equations in self.kernels_and_equations:
            if isinstance(kernel_or_equations, MatchedKernel):
                code_list.append(self._matched_kernel_to_generation_steps(kernel_or_equations))
            else:
                code_list.append(str(kernel_or_equations))

        code_list.append(str(self.final_equations))

        return "\n\n".join(code_list)

    def generation_steps_website(self):
        """Generates a description of how the algorithm was found for the website.

        Both equations and matched kernels are represented as strings containing
        tex code.

        Returns:
            list: A list of tuples of an integer and a string. The integer
                denotes whether the string contains an equation (0) or matched
                kernels (1).
        """
        output = []

        for is_matched_kernel, group in itertools.groupby(self.kernels_and_equations, lambda x: isinstance(x, MatchedKernel)):
            if is_matched_kernel:
                strings = []
                for matched_kernel in group:
                    strings.append(matched_kernel.to_tex(algin_char=True))
                tex_code = "\\begin{{align}}\n{}\n\\end{{align}}".format("\\\\\n".join(strings))
                output.append((1, tex_code)) 
            else:
                for equations in group:
                    output.append((0, equations.to_tex()))
        
        output.append((0, self.final_equations.to_tex()))
        return output

    def remove_duplicate_lines(self, lines, known_lines):
        """Removes all strings from lines that are in known_lines.

        Args:
            lines (list): List of strings.
            known_lines (set): Set of strings.

        Returns:
            list: The input list with known_lines removed.
        
        """
        new_lines = []
        for line in lines:
            if line not in known_lines:
                known_lines.add(line)
                new_lines.append(line)
        return new_lines

    def _generate_docstring(self, memory):
        input_operands, output_operands = self.initial_equations.input_output()
        operand_description = []
        internal_properties = {Property.AUXILIARY, Property.SCALAR,
                               Property.VECTOR, Property.MATRIX,
                               Property.SQUARE, Property.ROW_PANEL,
                               Property.COLUMN_PANEL, Property.CONSTANT,
                               Property.ADMITS_FACTORIZATION,
                               Property.FACTOR, Property.TRIANGULAR}
        for operand in input_operands:
            print_properties = operand.properties - internal_properties
            remove = set()
            # removing properties that are implied by others
            for p1, p2 in itertools.combinations(print_properties, 2):
                if p1 < p2:
                    remove.add(p2)
                elif p1 > p2:
                    remove.add(p1)
            print_properties -= remove
            if print_properties:
                properties_list = ", ".join(sorted([p.value for p in print_properties]))
                if len(print_properties) > 1:
                    properties_str = " with properties {}".format(properties_list)
                else:
                    properties_str = " with property {}".format(properties_list)
            else:
                properties_str = ""
            
            if isinstance(operand, Matrix):
                op_str = "Matrix {} of size {} x {}".format(operand.name, operand.rows, operand.columns)
            elif isinstance(operand, Vector):
                length = operand.columns if operand.columns != 1 else operand.rows
                op_str = "Vector {} of size {}".format(operand.name, length)
            elif isinstance(operand, Scalar):
                op_str = "Scalar {}".format(operand.name)

            description = "- `{}::{}`: {}{}.".format(memory.lookup[operand.name].name, operand_type(operand), op_str, properties_str)
            operand_description.append(description)
        template = textwrap.dedent(
                    """\
                    \"\"\"
                        {}({})

                    Compute
                    {}.

                    Requires at least Julia v1.0.

                    # Arguments
                    {}
                    \"\"\"\
                    """)
        return template.format(self.name, self._experiment_input, self.initial_equations, "\n".join(operand_description))

class KernelIO():
    """Describes the input and output of a kernel.

    This object describes the input and output of a kernel.

    Some remarks regarding the semantics of this object: For each argument (in
    the signature), it stores all the operands that are passed via this
    argument. Thus, one argument does not necessarily coincide with one operand.
    For factorizations, this is frequently the case, for example in LU, when
    both L and U are stored in the same array. Thus, they are returned by the
    function via the same argument.

    Consequently, in general, there is also no one-to-one correspondence to the
    substitution, because of functions to solve linear systems that expect for
    example one array containing both L and U as input. In that case, the
    substitution will contain two entries {L: L1, U: U2}, but the function only
    takes one argument for both (for example called A).

    In practice, for BLAS kernels, there is a one-to-one correspondence to the
    substitution, because one operand is always passed via one argument. Should
    we ever use kernels that work differently, then an additional step is
    necessary to get from the substitution to the input for add_input. It is
    questionable whether this will be supported anytime soon because it only
    makes sense to use such kernels if it is known that the operands are stored
    in the same memory location, which is currently not considered for pattern
    matching.

    Since those kernels are not supported yet, there is no list of operand for
    the input, only for the output.

    Attributes:
        entries (dict): Maps variables (that stand for arguments in the
            signature) to lists. Each list contains two elements. The first one
            is either a tuple (Expression, StorageFormat) of an input operand
            and its storage format, or None. The second element is a list of
            tuples (Expression, StorageFormat) of output operands and their
            storage formats.
        input_operands (list): Contains tuples (Expression, StorageFormat) of
            all input operands and their storage formats.
    """
    def __init__(self):
        self.entries = dict()
        self.input_operands = []

    def add_input(self, variable, operand, storage_format):
        input = (operand, storage_format)

        entry = self.entries.setdefault(variable, [None, []])
        if entry[0]:
            raise ValueError("More than one operand per input argument is currently not supported.")
        else:
            entry[0] = input

        self.input_operands.append(input)

    def add_output(self, variable, operand, storage_format):
        entry = self.entries.setdefault(variable, [None, []])
        entry[1].append((operand, storage_format))

    # def input_variables(self):
    #     variables = set()
    #     for variable, (input, output) in self.entries.items():
    #         if input:
    #             variables.add(variable.name)
    #     return variables

    def output_operands(self):
        for variable, (input, output) in self.entries.items():
            if input:
                overwritten_operand = input[0]
                for output_operand, storage_format in output:
                    yield (output_operand, overwritten_operand, storage_format)
            else:
                for output_operand, storage_format in output:
                    yield (output_operand, None, storage_format)

    def replace_variables(self, mapping):
        self.entries = {mapping[variable.name]: value for variable, value in self.entries.items()}

    def __repr__(self):
        str_list = []
        for variable, (input, output) in self.entries.items():
            entry_str = str(variable) + ","
            if input:
                entry_str = "".join([entry_str, " in: ", "({0}, {1})".format(str(input[0]), input[1].name)])
            if output:
                entry_str = "".join([entry_str, " out: ", self._io_to_str(output)])
            entry_str = "".join(["(", entry_str, ")"])
            str_list.append(entry_str)
        return "".join(["KernelIO(", ", ".join(str_list), ")"])

    def _io_to_str(self, io):
        str_list = []
        for operand, storage_format in io:
            str_list.append("({0}, {1})".format(str(operand), storage_format.name))
        return "".join(["(", ", ".join(str_list), ")"])

    # def __str__(self):
        # return str(self.entries)

    def __deepcopy__(self, memo):
        cpy = object.__new__(type(self))
        cpy.entries = copy.copy(self.entries)
        cpy.input_operands = self.input_operands.copy()
        return cpy


class MatchedKernel():
    """docstring for MatchedKernel"""

    _counter = 0

    def __init__(self):
        self.id = MatchedKernel._counter
        MatchedKernel._counter += 1
        
        self.operand_dict = None

        self.replacement = None
        self.operation = None
        self.cost = 0 # TODO change this to None. Currently, it's necessary to make sums work.

        self.signature = None
        self.arguments = None
        self.pre_code = None
        self.post_code = None

        self.other_replacements = None

        self.kernel_io = None

    def __eq__(self, other):
        return self.signature == other.signature and self.operation == other.operation

    def __hash__(self):
        return hash((self.signature, self.operation))

    def to_tex(self, algin_char=False):
        if algin_char:
            template = "{} &\\leftarrow {}"
        else:
            template = "{} \\leftarrow {}"
        return template.format(self.operation.lhs.to_tex(), self.operation.rhs.to_tex())


def to_file(file_path, content):
    directory_name = os.path.dirname(file_path)
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)
    output_file = open(file_path, "wt")
    output_file.write(content)
    output_file.close()
    if config.verbosity >= 2:
        print("Generate file {}".format(file_path))

def operand_type(operand, property_types=False):
    if property_types:
        if operand.has_property(Property.SCALAR):
            return config.data_type_string
        elif operand.has_property(Property.VECTOR):
            return "Array{{{0},1}}".format(config.data_type_string)
        elif operand.has_property(Property.DIAGONAL):
            return "Diagonal{{{0},Array{{{0},1}}}}".format(config.data_type_string)
        elif operand.has_property(Property.SYMMETRIC) or operand.has_property(Property.SPD) or operand.has_property(Property.SPSD):
            return "Symmetric{{{0},Array{{{0},2}}}}".format(config.data_type_string)
        elif operand.has_property(Property.LOWER_TRIANGULAR) and operand.has_property(Property.SQUARE):
            return "LowerTriangular{{{0},Array{{{0},2}}}}".format(config.data_type_string)
        elif operand.has_property(Property.UPPER_TRIANGULAR) and operand.has_property(Property.SQUARE):
            return "UpperTriangular{{{0},Array{{{0},2}}}}".format(config.data_type_string)
        else:
            return "Array{{{0},2}}".format(config.data_type_string)
    else:
        if operand.has_property(Property.SCALAR):
            return config.data_type_string
        elif operand.has_property(Property.VECTOR):
            return "Array{{{0},1}}".format(config.data_type_string)
        else:
            return "Array{{{0},2}}".format(config.data_type_string)

def remove_files(directory_name):
    if os.path.exists(directory_name):
        for file in os.listdir(directory_name):
            path_to_file = os.path.join(directory_name, file)
            if os.path.isfile(path_to_file):
                os.remove(path_to_file)


WD1 = matchpy.Wildcard.dot("WD1")
WD2 = matchpy.Wildcard.dot("WD2")
WS1 = matchpy.Wildcard.star("WS1")
WS2 = matchpy.Wildcard.star("WS2")
PS1 = matchpy.CustomConstraint(lambda WD1: WD1.has_property(Property.MATRIX) or WD1.has_property(Property.VECTOR))
PS2 = matchpy.CustomConstraint(lambda WD2: WD2.has_property(Property.MATRIX) or WD2.has_property(Property.VECTOR))
notInv1 = matchpy.CustomConstraint(lambda WD1: not is_inverse(WD1))
notInv2 = matchpy.CustomConstraint(lambda WD2: not is_inverse(WD2))

linsolveL = matchpy.ReplacementRule(
    matchpy.Pattern(Times(WS1, Inverse(WD1), WD2, WS2), PS1, PS2),
    lambda WS1, WD1, WD2, WS2: Times(*WS1, LinSolveL(WD1, WD2), *WS2)
    )

linsolveLT = matchpy.ReplacementRule(
    matchpy.Pattern(Times(WS1, InverseTranspose(WD1), WD2, WS2), PS1, PS2),
    lambda WS1, WD1, WD2, WS2: Times(*WS1, LinSolveL(Transpose(WD1), WD2), *WS2)
    )

linsolveR = matchpy.ReplacementRule(
    matchpy.Pattern(Times(WS1, WD1, Inverse(WD2), WS2), PS1, PS2),
    lambda WS1, WD1, WD2, WS2: Times(*WS1, LinSolveR(WD1, WD2), *WS2)
    )

linsolveRT = matchpy.ReplacementRule(
    matchpy.Pattern(Times(WS1, WD1, InverseTranspose(WD2), WS2), PS1, PS2),
    lambda WS1, WD1, WD2, WS2: Times(*WS1, LinSolveR(WD1, Transpose(WD2)), *WS2)
    )

linsolveLnI = matchpy.ReplacementRule(
    matchpy.Pattern(Times(WS1, Inverse(WD1), WD2, WS2), PS1, PS2, notInv2),
    lambda WS1, WD1, WD2, WS2: Times(*WS1, LinSolveL(WD1, WD2), *WS2)
    )

linsolveLTnI = matchpy.ReplacementRule(
    matchpy.Pattern(Times(WS1, InverseTranspose(WD1), WD2, WS2), PS1, PS2, notInv2),
    lambda WS1, WD1, WD2, WS2: Times(*WS1, LinSolveL(Transpose(WD1), WD2), *WS2)
    )

linsolveRnI = matchpy.ReplacementRule(
    matchpy.Pattern(Times(WS1, WD1, Inverse(WD2), WS2), PS1, PS2, notInv1),
    lambda WS1, WD1, WD2, WS2: Times(*WS1, LinSolveR(WD1, WD2), *WS2)
    )

linsolveRTnI = matchpy.ReplacementRule(
    matchpy.Pattern(Times(WS1, WD1, InverseTranspose(WD2), WS2), PS1, PS2, notInv1),
    lambda WS1, WD1, WD2, WS2: Times(*WS1, LinSolveR(WD1, Transpose(WD2)), *WS2)
    )

def replace_linsolve_left(equations):
    """Replaces linear systems with a LinSolve operator.

    This function can be used to generate Equations to generate naive and
    recommended code.
    """
    new_eqns = []
    for eqn in equations:
        eqn = matchpy.replace_all(eqn, [linsolveLnI, linsolveLTnI])
        eqn = matchpy.replace_all(eqn, [linsolveL, linsolveLT])
        new_eqns.append(eqn)
    return Equations(*new_eqns)

    
def replace_linsolve_right(equations):
    """Replaces linear systems with a LinSolve operator.

    This function can be used to generate Equations to generate naive and
    recommended code.
    """
    new_eqns = []
    for eqn in equations:
        eqn = matchpy.replace_all(eqn, [linsolveRnI, linsolveRTnI])
        eqn = matchpy.replace_all(eqn, [linsolveR, linsolveRT])
        new_eqns.append(eqn)
    return Equations(*new_eqns)


def get_template(file_name, language):
    template_path = "experiments/templates/{}/{}".format(language.name, file_name)
    return pkg_resources.resource_string(__name__, template_path).decode("UTF-8")