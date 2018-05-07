from ...algebra.expression import Matrix, Equal, Times
from ...algebra.transformations import transpose, conjugate_transpose, invert
from ...algebra.properties import Property as properties

from ... import temporaries

from ..utils.general import substitute_symbols_with_wildcards, \
                               to_wildcard_name, \
                               to_c_variable, to_c_variable_definition, \
                               SizeArgument, \
                               Kernel, \
                               InputOperand

from ...code_generation.utils import MatchedKernel, KernelIO

from ... import config

import copy
import operator
import itertools

import matchpy

class FactorizationKernel(Kernel):
    """docstring for FactorizationKernel"""

    def __init__(self, pattern, input_operands, output, output_operands, cost_function, pre_code, signature, post_code, arguments):
        super().__init__(cost_function, pre_code, signature, post_code, arguments)

        self.pattern = pattern
        self.replacement_template = output
        self.input_operands = input_operands
        self.output_operands = output_operands
        self.operation_template = self.pattern

        self.pre_code_template = pre_code
        self.post_code_template = post_code

        # TODO for something like generalized schur decomposition, I potentially also need context

    def set_match(self, match_dict, context, CSE_rules=False, set_equivalent=True, equiv_expr=None):

        matched_kernel = super().set_match(match_dict, CSE_rules)

        #############
        # operation

        # Constructing the input expression
        _input_expr = matchpy.substitute(self.operation_template, match_dict)

        # create output expression and arg_dict
        try:
            op_dict = temporaries._table_of_factors[self.id]
        except KeyError:
            # if there is no dict for current factorization, create
            # everything and store them in new dict
            ops = self._set_match(match_dict, _input_expr)
            temporaries._table_of_factors[self.id] = {_input_expr: ops}
        else:
            try:
                ops = op_dict[_input_expr]
            except KeyError:
                # if there is nothing stored for this match, create and store everything
                ops = self._set_match(match_dict, _input_expr)
                op_dict[_input_expr] = ops

        _output_expr, _arg_dict, _partial_operand_dict, kernel_io = ops
        # print(_partial_operand_dict)

        matched_kernel.operation = Equal(_output_expr, _input_expr)

        #############
        # operand_dict

        matched_kernel.kernel_io = kernel_io

        matched_kernel.operand_dict = copy.copy(match_dict)
        matched_kernel.operand_dict.update(_partial_operand_dict)
        # print(matched_kernel.operand_dict)
        # print(match_dict, _output_expr, _output, matched_kernel.operand_dict)

        #############
        # Replacement

        if context:
            # When this gets implemented, don't forget to remove context variables from match_dict
            raise NotImplementedError()
        else:
            _replacement = _output_expr

        matched_kernel.replacement = _replacement

        #############
        # Other replacements

        matched_kernel.other_replacements = {
                    "type": config.data_type_string,
                    "type_prefix": config.blas_data_type_prefix, # TODO this is language dependent
                    "work_id": hex(hash(self.signature))[-5:]
                    }

        #############
        # CSE rules

        if CSE_rules:
            _rules = [
                (
                    matchpy.Pattern(_input_expr),
                    lambda **_: _output_expr
                ),
                (
                    matchpy.Pattern(transpose(_input_expr)),
                    lambda **_: transpose(_output_expr)
                ),
                (
                    matchpy.Pattern(conjugate_transpose(_input_expr)),
                    lambda **_: conjugate_transpose(_output_expr)
                ),
            ]
            matched_kernel.CSE_rules = _rules

        return matched_kernel

    def _set_match(self, match_dict, input_expr):
        """Auxiliary function for set_match()

        Computes only those things that are independent of whether temporaries
        are reused or not.
        """

        # Constructing input.
        kernel_io = KernelIO()
        for input_operand in self.input_operands:
            kernel_io.add_input(input_operand.operand, match_dict[input_operand.operand.variable_name], input_operand.storage_format)

        _arg_dict = dict()
        for arg in self.arguments:
            if isinstance(arg, SizeArgument):
                _arg_dict[arg.name] = arg.get_value(match_dict)

        _partial_operand_dict = dict()

        # replacement_dict maps wildcard names to operands
        replacement_dict = dict()
        for output_operand in self.output_operands:
            # output_operand.operand.name[1:] because it's a Wildcard, we drop the _
            name = "".join([output_operand.operand.variable_name[1:], temporaries.get_identifier()])
            size = (_arg_dict[output_operand.size[0]], _arg_dict[output_operand.size[1]])
            # TODO what if the output is a scalar? Check sizes.

            operand = Matrix(name, size, input_expr.indices)
            operand.set_property(properties.FACTOR)
            operand.factorization_labels = set(operand[0].name for operand in kernel_io.input_operands)
            for property in output_operand.properties:
                operand.set_property(property)

            replacement_dict[output_operand.operand.variable_name] = operand


            # Constructing output.
            if output_operand.overwriting:
                kernel_io.add_output(output_operand.overwriting, operand, output_operand.storage_format)
            else:
                kernel_io.add_output(output_operand.operand, operand, output_operand.storage_format)
                _partial_operand_dict[output_operand.operand.variable_name] = operand

        _output_expr = matchpy.substitute(self.replacement_template, replacement_dict)

        return _output_expr, _arg_dict, _partial_operand_dict, kernel_io


class OutputOperand():
    """docstring for OutputOperand"""
    def __init__(self, operand, overwriting, size, properties, storage_format):
        self.operand = operand
        self.overwriting = overwriting
        self.size = size
        self.properties = properties
        self.storage_format = storage_format
        # TODO missing: what is the corresponding argument in the function signature? Only necessary of not overwriting.

    def __repr__(self):
        return "".join(["OutputOperand(",
                            self.operand.name, ", ",
                            self.overwriting.name, ", ",
                            self.size, ", ",
                            self.properties, ", ",
                            self.storage_format.name, ")"])

