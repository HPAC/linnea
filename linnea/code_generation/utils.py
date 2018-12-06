
from ..algebra.expression import Symbol, Times, \
                                 Inverse, InverseTranspose, Transpose, \
                                 LinSolveL, LinSolveR
from ..algebra.equations import Equations

from .. import config

from ..utils import PropertyConstraint

from .memory import memory as memory_module

from ..algebra.properties import Property as properties

from ..derivation.graph.utils import is_inverse

import copy
import math
import textwrap
import os
import matchpy
import pkg_resources

class Algorithm():
    """Represents an Algorithm and translates it to code.

    This class represents one algorithm and offers functions to translate it to
    code.

    Attributes:
        initial_equations (Equations): The input equations.
        final_equations (Equations): The equations left at the end of the
            dervation. They only consist of quations of the form
            Symbol = Symbol. Is used for the mapping of temporaries to output
            operands.
        kernels_and_equations (list): List of MatchedKernel and Equation
            objects, representing the algorithm (and its derivation).
        cost (int): The cost of the entire algorithm.
        memory (Memory): A Memory object, which is used to translate the
            algorithm to code.
        liveness (dict): Stores liveness information. For each operand, it
            contains a tuple of to integers. The first integer is the line
            number (corresponds to a MatchedKernel object in matched_kernels)
            where the operand is define. If it is None, the operand is an input
            operand. The second integer is the line number where the operand is
            used for the last time. If it is None, the operand is an ouput
            operand.
        known_lines (set): Is used to generate code. Contains lines of code that
            may be generated multiple times even though they are only needed
            once, such as definitions of constant arguments.
    """
    def __init__(self, initial_equations, final_equations, kernels_and_equations, cost):

        # super(Algorithm, self).__init__()
        self.initial_equations = initial_equations
        self.final_equations = final_equations
        self.matched_kernels = [keq for keq in kernels_and_equations if isinstance(keq, MatchedKernel)]
        self.kernels_and_equations = kernels_and_equations
        self.cost = cost

        self.memory = None

        self.liveness = None
        self.known_lines = None

        self.experiment_input = None
        self.experiment_output = None

        self.liveness_analysis()

    def _symbols(self, expr):
        """Returns all the symbols in expr.
    
        """
        if isinstance(expr, Symbol):
            yield expr
        else:
            for operand in expr.operands:    
                yield from self._symbols(operand)

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
        self.liveness = dict()
        live_symbols = set()
        for line_number, matched_kernel in reversed(list(enumerate(self.matched_kernels))):
            rhs_symbols = set(self._symbols(matched_kernel.operation.rhs))
            for rhs_symbol in rhs_symbols:
                if rhs_symbol not in live_symbols:
                    self.liveness.setdefault(rhs_symbol.name ,[None, None])[1] = line_number
            
            live_symbols.update(rhs_symbols)

            lhs_symbols = set(self._symbols(matched_kernel.operation.lhs))
            for lhs_symbol in lhs_symbols:
                self.liveness.setdefault(lhs_symbol.name ,[None, None])[0] = line_number

        for equation in self.final_equations:
            self.liveness[equation.rhs.name][1] = None

    def code(self):
        """Translated the algorithm to code.

        The language depends on what is specified in config.
        
        Returns:
            string: Code.
        """


        """
        TODO this is inconsistent. memory and known_lines are attributes,
        matched_kernel and line_number is passed to memory.add_operation
        """

        self.memory = memory_module.Memory(self.initial_equations)
        self.known_lines = set()
        code_list = []

        input_operands, output_operands = self.initial_equations.input_output()
        self.experiment_input = ", ".join([self.memory.lookup[operand.name].name for operand in input_operands])

        code_list.append("{0}cost {1:.3g}\n".format(config.comment, self.cost))

        if config.c:
            code_list.append("int info = 0;\n\n")

        for line_number, matched_kernel in enumerate(self.matched_kernels):
            code_list.append(self._matched_kernel_to_code(matched_kernel, line_number))

        code_list.append(config.comment)
        code_list.append(self.memory.content_string_with_format())
        code_list.append("\n")

        output_operand_mapping = {eqn.lhs: eqn.rhs for eqn in self.final_equations}
        conversion_to_full = self.memory.convert_to_full(output_operand_mapping.values())
        if conversion_to_full:
            code_list.append(conversion_to_full)
            code_list.append("\n")

            code_list.append(config.comment)
            code_list.append(self.memory.content_string_with_format())
            code_list.append("\n")

        code_list.append(textwrap.indent(str(self.final_equations), config.comment))
        self.experiment_output = ", ".join([self.memory.lookup[output_operand_mapping[operand].name].name for operand in output_operands])

        return "".join(code_list)

    def _matched_kernel_to_code(self, matched_kernel, line_number):
        """Generates code for a single MatchedKernel objects.

        Args:
            matched_kernel (MatchedKernel)
            line_number (int): The line number of the operation represeted by
                matched_kernel. Is necessary to use the liveness dictionary.

        Returns:
            string: Code.

        Important:
            MatchedKernel objects can not be modified because they are used
            multiple times for the generation of algorithms.

        """
        lines_list = []

        # print("#########")
        # print(matched_kernel.signature)
        # print(matched_kernel.operand_dict)
        # print(matched_kernel.kernel_io)

        mem_content = "".join([config.comment, self.memory.content_string_with_format(), "\n"])
        # mem_content = "".join([config.comment, self.memory.content_string(), "\n"])
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
            try:
                pre_code, arg_replacement = argument.get_replacement(matched_kernel.operand_dict, self.memory)
            except memory_module.OperandNotInMemory:
                """ This Argument object has to be used after
                self.memory.add_operation(). Happens only for StrideArgument objects
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
        mem_ops_before, mem_ops_after, operand_mapping = self.memory.add_operation(matched_kernel.kernel_io, self.liveness, line_number)
        # print(operand_mapping)

        if mem_ops_before:
            mem_code_before = "".join([mem_op.code() for mem_op in mem_ops_before])
            lines_list.append(mem_code_before)

        # print(late_arguments)
        for argument in late_arguments:
            pre_code, arg_replacement = argument.get_replacement(matched_kernel.operand_dict, self.memory)
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
            lines_list.extend(self.remove_duplicate_lines(argument_pre_code, self.known_lines))

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


    def _matched_kernel_to_derivation(self, matched_kernel):
        """Generates pseudocode for a signle MatchedKernel objects.

        Args:
            matched_kernel (MatchedKernel)

        Important:
            MatchedKernel objects can not be modified because they are used
            multiple times for the generation of algorithms.

        Returns:
            string: Pseudocode.
        """
        return "{0:<30}# {1:.3g}".format(str(matched_kernel.operation), matched_kernel.cost)

    def derivation(self):
        """Generates a description of how the algorithm was found.
        
        Returns:
            string: Derivation.
        """
        code_list = []

        code_list.append("# cost {:.3g}".format(self.cost))

        for kernel_or_equations in self.kernels_and_equations:
            if isinstance(kernel_or_equations, MatchedKernel):
                code_list.append(self._matched_kernel_to_derivation(kernel_or_equations))
            else:
                code_list.append(str(kernel_or_equations))

        code_list.append(str(self.final_equations))

        return "\n\n".join(code_list)

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

    def __init__(self, CSE_rules=False):
        self.id = MatchedKernel._counter
        MatchedKernel._counter += 1
        
        self.operand_dict = None

        self.replacement = None
        self.operation = None
        self.cost = 0 # TODO change this to None. Currently, it's necessary to make sums work.

        if CSE_rules:
            self.CSE_rules = None

        self.signature = None
        self.arguments = None
        self.pre_code = None
        self.post_code = None

        self.other_replacements = None

        self.kernel_io = None
        

def algorithm_to_file(output_name, subdir_name, algorithm_name, algorithm, input, output,
                      language = config.language,
                      file_extension = config.filename_extension):
    file_name = os.path.join(config.output_code_path, output_name, language.name, subdir_name,
                             "".join([algorithm_name, file_extension]))
    directory_name = os.path.dirname(file_name)
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)
    output_file = open(file_name, "wt")
    if config.verbosity >= 2:
        print("Generate algorithm file {}".format(file_name))
    algorithm_name = algorithm_name
    algorithm_str = algorithm_to_str(algorithm_name, algorithm, input, output, language)
    output_file.write(algorithm_str)
    output_file.close()

def algorithm_to_str(function_name, algorithm, input, output, language):
    if language == config.Language.Julia:
        template = get_template("algorithm.jl", language)
        algorithm_str = template.format(function_name, input, textwrap.indent(algorithm, "    "), output)
    elif language == config.Language.Matlab:
        template = get_template("algorithm.m", language)
        algorithm_str = template.format(function_name, input, textwrap.indent(algorithm, "    "), output)
    elif language == config.Language.Cpp:
        types_list = ", ".join("typename Type_{}".format(op) for op in str.split(input, ", "))
        args_list = ", ".join("Type_{0} && {0}".format(op) for op in str.split(input, ", "))
        output_len = len(str.split(output, ", "))
        if output_len > 1:
            output_string = "return std::make_tuple({});".format(output)
        else:
            ret_string = textwrap.dedent(
                         """
                            typedef std::remove_reference_t<decltype({})> return_t;
                            return return_t({});
                         """)
            output_string = ret_string.format(output, output)
        template = get_template("algorithm.cpp", language)
        algorithm_str = template.format(function_name, types_list, args_list, textwrap.indent(algorithm, "    "), output_string)
    else:
        raise config.LanguageOptionNotImplemented()
    
    return algorithm_str

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
PS1 = matchpy.CustomConstraint(lambda WD1: WD1.has_property(properties.MATRIX) or WD1.has_property(properties.VECTOR))
PS2 = matchpy.CustomConstraint(lambda WD2: WD2.has_property(properties.MATRIX) or WD2.has_property(properties.VECTOR))
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