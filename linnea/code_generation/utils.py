
from ..algebra.expression import Symbol

from .. import config

from .memory import memory as memory_module

import copy
import math
import textwrap

class Algorithm(object):
    """Represents an Algorithm and translates it to code.

    This class represents one algorithm and offers functions to translate it to
    code.

    Attributes:
        initial_equations (Equations): The input equations.
        final_equations (Equations): The equations left at the end of the
            dervation. They only consist of quations of the form
            Symbol = Symbol. Is used for the mapping of temporaries to output
            operands.
        matched_kernels (list): List of MatchedKernel objects, representing the
            operations of the algorithm.
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
    def __init__(self, initial_equations, final_equations, matched_kernels, cost):

        # super(Algorithm, self).__init__()
        self.initial_equations = initial_equations
        self.final_equations = final_equations
        self.matched_kernels = matched_kernels
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
        elif config.julia:
            code_list.append("using Base.LinAlg.BLAS\nusing Base.LinAlg\n\n")

        for line_number, matched_kernel in enumerate(self.matched_kernels):
            code_list.append(self._matched_kernel_to_code(matched_kernel, line_number))

        code_list.append(config.comment)
        code_list.append(self.memory.content_string_with_format())
        code_list.append("\n")

        code_list.append(textwrap.indent(str(self.final_equations), config.comment))

        output_operand_mapping = {eqn.lhs: eqn.rhs for eqn in self.final_equations}
        self.experiment_output = ", ".join([self.memory.lookup[output_operand_mapping[operand].name].name for operand in output_operands])

        return "".join(code_list)

    def _matched_kernel_to_code(self, matched_kernel, line_number):
        """Generates code for a signle MatchedKernel objects.

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
        lines = ""

        # print("#########")
        # print(matched_kernel.signature)
        # print(matched_kernel.operand_dict)
        # print(matched_kernel.kernel_io)

        mem_content = "".join([config.comment, self.memory.content_string_with_format(), "\n"])
        # mem_content = "".join([config.comment, self.memory.content_string(), "\n"])
        lines = "".join([lines, mem_content])

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
            lines = "".join([lines, mem_code_before])

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
            lines = "".join([lines, *matched_kernel.remove_duplicate_lines(argument_pre_code, self.known_lines)])

        if matched_kernel.pre_code:
            lines = "".join([lines, kernel_pre_code.safe_substitute_str(argument_mapping)])

        lines = "".join([lines, config.comment, str(matched_kernel.operation), "\n"])
        lines = "".join([lines, signature.safe_substitute_str(argument_mapping), "\n"])

        if matched_kernel.post_code:
            lines = "".join([lines, kernel_post_code.safe_substitute_str(argument_mapping)])

        if mem_ops_after:
            mem_code_after = "".join([mem_op.code() for mem_op in mem_ops_after])
            lines = "".join([lines, mem_code_after])

        lines = "".join([lines, "\n"])
        return lines


    def _matched_kernel_to_pseudocode(self, matched_kernel):
        """Generates pseudocode for a signle MatchedKernel objects.

        Args:
            matched_kernel (MatchedKernel)

        Important:
            MatchedKernel objects can not be modified because they are used
            multiple times for the generation of algorithms.

        Returns:
            string: Pseudocode.
        """
        return "{0:<30}# {1:.3g}\n".format(str(matched_kernel.operation), matched_kernel.cost)

    def pseudocode(self):
        """Translated the algorithm to pseudocode.
        
        Returns:
            string: Pseudocode.
        """
        code_list = []

        code_list.append("# cost {:.3g}\n".format(self.cost))

        for matched_kernel in self.matched_kernels:
            code_list.append(self._matched_kernel_to_pseudocode(matched_kernel))

        code_list.append(str(self.final_equations))

        return "".join(code_list)

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

class KernelIO(object):
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


class MatchedKernel(object):
    """docstring for MatchedKernel"""
    def __init__(self, CSE_rules=False):
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
        
    def code(self, memory, known_lines):
        """
        
        Important:
            This function can not modify self. The reason is that it may get
            called multiple times, for different algorithms, with different
            input arguments.
        """
        lines = ""

        # print("#########")
        # print(self.signature)
        # print(self.operand_dict)
        # print(self.kernel_io)

        mem_content = "".join([config.comment, memory.content_string(), "\n"])
        lines = "".join([lines, mem_content])

        # TODO arguments could be a stack or something, removing processed arguments. Then we don't need late_arguments.
        arguments = copy.copy(self.arguments)
        argument_mapping = dict()
        argument_pre_code = []

        signature = self.signature.safe_substitute_copy(self.other_replacements)
        if self.pre_code:
            self.pre_code.safe_substitute(self.other_replacements)
        if self.post_code:
            self.post_code.safe_substitute(self.other_replacements)

        late_arguments = []
        for argument in arguments:
            try:
                pre_code, arg_replacement = argument.get_replacement(self.operand_dict, memory)
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

        mem_ops_before, mem_ops_after, operand_mapping = memory.add_operation(self.kernel_io)
        # print(operand_mapping)

        if mem_ops_before:
            mem_code_before = "".join([mem_op.code() for mem_op in mem_ops_before])
            lines = "".join([lines, mem_code_before])

        # print(late_arguments)
        for argument in late_arguments:
            pre_code, arg_replacement = argument.get_replacement(self.operand_dict, memory)
            if pre_code:
                argument_pre_code.append(pre_code)
            argument_mapping[argument.name] = arg_replacement 

        # print(argument_mapping)
        # print(signature)

        signature.safe_substitute(operand_mapping)
        if self.pre_code:
            self.pre_code.safe_substitute_copy(operand_mapping)
        if self.post_code:
            self.post_code.safe_substitute_copy(operand_mapping)

        if argument_pre_code:
            lines = "".join([lines, *self.remove_duplicate_lines(argument_pre_code, known_lines)])

        if self.pre_code:
            lines = "".join([lines, self.pre_code.safe_substitute_str(argument_mapping)])

        lines = "".join([lines, config.comment, str(self.operation), "\n"])
        lines = "".join([lines, signature.safe_substitute_str(argument_mapping), "\n"])

        if self.post_code:
            lines = "".join([lines, self.post_code.safe_substitute_str(argument_mapping)])

        if mem_ops_after:
            mem_code_after = "".join([mem_op.code() for mem_op in mem_ops_after])
            lines = "".join([lines, mem_code_after])

        lines = "".join([lines, "\n"])
        return lines


    def pseudocode(self):
        return "{0:<30}# {1:.0f}\n".format(str(self.operation), self.cost)

    def remove_duplicate_lines(self, lines, known_lines):
        new_lines = []
        for line in lines:
            if line not in known_lines:
                known_lines.add(line)
                new_lines.append(line)
        return new_lines

