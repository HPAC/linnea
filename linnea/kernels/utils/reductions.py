from ...algebra.expression import Equal, Times, Plus, \
                                  Operator, Symbol, \
                                  Inverse, Transpose, InverseTranspose

from ...algebra.transformations import simplify, \
                                       transpose, conjugate_transpose

from .general import PropertyConstraints, \
                     to_c_variable, to_c_variable_definition, \
                     to_wildcard_name, \
                     SizeArgument, ConstantArgument, \
                     substitute_symbols_with_wildcards, \
                     Kernel, \
                     InputOperand, \
                     PropertyTuple

from ...utils import CodeTemplate, PropertyConstraint

from ... import config

from ...code_generation.memory import storage_format
from ...code_generation.utils import MatchedKernel, KernelIO

from ...algebra.properties import Property as properties

from ... import temporaries

import enum
import copy
import itertools
import string
import operator
import matchpy

class KernelType(enum.Enum):
    identity = 0
    transpose = 1
    conjugate_transpose = 2

_op = matchpy.Wildcard.symbol("_op")
ctx1 = matchpy.Wildcard.star("ctx1")
ctx2 = matchpy.Wildcard.star("ctx2")

class ReductionKernel(Kernel):
    """docstring for ReductionKernel"""

    def __init__(self, pattern, input_operands, output_operand, cost_function, pre_code, signature, post_code, arguments, type=KernelType.identity, kernel_io=None):
        super().__init__(cost_function, pre_code, signature, post_code, arguments)

        """The wildcard which represents the output operand. If it also appears
        in the pattern, it is assumed that it gets overwritten.
        """
        self.output_operand = output_operand

        self.input_operands = input_operands

        self.is_overwriting = False
        for input_operand in self.input_operands:
            if input_operand.operand == self.output_operand.operand:
                self.is_overwriting = True
                break

        self.kernel_io = kernel_io
        if not self.kernel_io:
            self.kernel_io = KernelIO()

        # print(self.input_operands)
        # print(self.output_operand)

        # Template for the operation that is computed. Is used to construct the
        # temporary.
        # Does not contain context.
        # Does not take type into account.
        # Example TRSV (transpose type):
        #   Times(Inverse(X_), y_)
        #   (X is wildcard for a matrix, y for a vector)
        self.operation_template = pattern

        # The pattern object to be used for matching.
        # Can be used for the matrix chain algorithm.
        # Does not contain context.
        # Takes type into account.
        # Example TRSV (transpose type):
        #   Times(Transpose(y_), InverseTranspose(X_))
        self.pattern = pattern

        # Template for the replacement.
        # Does not contain context.
        # Takes type into account.
        # Example TRSV (transpose type):
        #   Transpose(op)
        self.replacement_template = _op

        # TODO is it necessary to simplify?

        pattern_expr = self.pattern.expression

        # Modifying pattern and replacement to take type into account.
        # TODO this doesn't work for some kernels (e.g. diagsv), because
        # variables don't have properties that might be necessary to simplify
        # expressions. Those properties are part of the constraints, not the
        # expression itself.
        if type == KernelType.identity:
            pass
        elif type == KernelType.transpose:
            pattern_expr = transpose(pattern_expr)
            # TODO ugly hack
            for constraint in self.pattern.constraints:
                if (properties.SQUARE in constraint.properties and properties.DIAGONAL in constraint.properties) or properties.SYMMETRIC in constraint.properties:
                    pattern_expr = remove_transpose(pattern_expr, constraint.variable)
            self.replacement_template = transpose(self.replacement_template)
        elif type == KernelType.conjugate_transpose:
            pattern_expr = conjugate_transpose(pattern_expr)
            self.replacement_template = conjugate_transpose(self.replacement_template)

        self.pattern = matchpy.Pattern(pattern_expr, *self.pattern.constraints)

        # The pattern object to be used for matching.
        # Contains context.
        # Takes type into account.
        # Example TRSV (transpose type):
        #   Times(ctx1___, Transpose(y_), InverseTranspose(X_), ctx2___)
        self.pattern_with_context = pattern

        # Template for the replacement.
        # Contains context.
        # Takes type into account.
        # Example TRSV (transpose type):
        #   Times(ctx1, Transpose(op), ctx2)
        self.replacement_with_context_template = None

        # adding contexts
        if isinstance(pattern.expression, Times):
            pattern_with_context_expr = Times(ctx1, self.pattern.expression, ctx2)
            self.replacement_with_context_template = Times(ctx1, self.replacement_template, ctx2)
        elif isinstance(pattern.expression, Plus):
            pattern_with_context_expr = Plus(ctx1, self.pattern.expression)
            self.replacement_with_context_template = Plus(ctx1, self.replacement_template)
        else: # For unary kernels
            pattern_with_context_expr = self.pattern.expression
            self.replacement_with_context_template = self.replacement_template

        self.pattern_with_context = matchpy.Pattern(pattern_with_context_expr, *self.pattern.constraints)

        property_list = []
        for expr, pos in self.pattern.expression.preorder_iter():
            if isinstance(expr, matchpy.Wildcard):
                for constraint in self.pattern.constraints:
                    if constraint.variable == expr.variable_name:
                        property_list.append(constraint.properties)
                        break
                else:
                    property_list.append(frozenset())
        self.property_tuple = PropertyTuple(property_list)
        # print(self.pattern, self.property_tuple)

        # print("here", self.id, self.pattern)

    # @profile
    def set_match(self, match_dict, context, CSE_rules=False, blocked_products=False, set_equivalent=True, equiv_expr=None):
        
        matched_kernel = super().set_match(match_dict, CSE_rules)

        #############
        # operation

        # I don't like this part. I would prefer not to use it at all and always
        # use equivalent expression (even for matrix chain). It's exclusively
        # a performance consideration.
        if equiv_expr:
            # equiv_expr = self.replacement_template.replace_copy({"_op": equiv_expr})    
            equiv_expr = matchpy.substitute(self.replacement_template, {"_op": equiv_expr})
            equiv_expr = simplify(equiv_expr)

        # _operation = self.operation_template.replace_copy(match_dict)
        _operation = matchpy.substitute(self.operation_template, match_dict)
        _tmp = temporaries.create_tmp(_operation, set_equivalent, equiv_expr)

        matched_kernel.operation = Equal(_tmp, _operation)

        #############
        # operand_dict & kernel_io

        kernel_io = copy.deepcopy(self.kernel_io)

        operand_dict = dict()
        for input_operand in self.input_operands:
            matched_operand = match_dict[input_operand.operand.variable_name]
            operand_dict[input_operand.operand.variable_name] = matched_operand
            kernel_io.add_input(input_operand.operand, matched_operand, input_operand.storage_format)

        kernel_io.add_output(self.output_operand.operand, _tmp, self.output_operand.storage_format)
        if not self.is_overwriting:
            operand_dict[self.output_operand.operand.variable_name] = _tmp

        matched_kernel.kernel_io = kernel_io

        matched_kernel.operand_dict = operand_dict


        #############
        # Replacement

        if context:
            # _replacement = matchpy.substitute(self.replacement_with_context_template, {"_op": _tmp}+match_dict)[0] # Plugging in tmp
            _replacement = matchpy.substitute(self.replacement_with_context_template, {"_op": _tmp}) # Plugging in tmp
            _replacement = matchpy.substitute(_replacement, match_dict) # Plugging in context
        else:
            # _replacement = self.replacement_template.replace_copy({"_op": _tmp})
            _replacement = matchpy.substitute(self.replacement_template, {"_op": _tmp})

        matched_kernel.replacement = _replacement

        #############
        # Other replacements

        # TODO This is language dependent
        matched_kernel.other_replacements = {
            "type": config.data_type_string,
            "type_prefix": config.blas_data_type_prefix,
            }

        #############
        # CSE rules

        if CSE_rules:
            if isinstance(self.pattern.expression, Times):
                _rules = [
                    (
                        matchpy.Pattern(Times(ctx1, _operation, ctx2)),
                        lambda ctx1, ctx2: Times(*ctx1, _tmp, *ctx2)
                    ),
                    (
                        matchpy.Pattern(Times(ctx1, transpose(_operation), ctx2)),
                        lambda ctx1, ctx2: Times(*ctx1, transpose(_tmp), *ctx2)
                    ),
                    (
                        matchpy.Pattern(Times(ctx1, conjugate_transpose(_operation), ctx2)),
                        lambda ctx1, ctx2: Times(*ctx1, conjugate_transpose(_tmp), *ctx2)
                    ),
                ]
            elif isinstance(self.pattern.expression, Plus):
                _rules = [
                    (
                        matchpy.Pattern(Plus(ctx1, _operation)),
                        lambda ctx1: Plus(*ctx1, _tmp)
                    ),
                    (
                        matchpy.Pattern(Plus(ctx1, transpose(_operation))),
                        lambda ctx1: Plus(*ctx1, transpose(_tmp))
                    ),
                    (
                        matchpy.Pattern(Plus(ctx1, conjugate_transpose(_operation))),
                        lambda ctx1: Plus(*ctx1, conjugate_transpose(_tmp))
                    ),
                ]
            else:
                _rules = [
                    (
                        matchpy.Pattern(_operation),
                        lambda **_: _tmp
                    ),
                    (
                        matchpy.Pattern(transpose(_operation)),
                        lambda **_: transpose(_tmp)
                    ),
                    (
                        matchpy.Pattern(conjugate_transpose(_operation)),
                        lambda **_: conjugate_transpose(_tmp)
                    ),
                ]
            matched_kernel.CSE_rules = _rules

        #############
        # Blocked products

        # Not relevant for reductions.   

        return matched_kernel

    def is_matrix_chain_kernel(self):
        """Returns True if the kernel can be used for the matrix chain algorithm."""
        pattern_expr = self.pattern.expression
        if isinstance(pattern_expr, Times) and len(pattern_expr.operands) == 2:
            return True
        else:
            return False

def remove_transpose(expr, name):
    """Removes transpose operator for a symbol.

    This is an ugly hack to temporarily solve the problem that simply doesn't
    work correctly for pattern because it can't consider properties.

    expr = remove_transpose(expr)
    """
    if isinstance(expr, Transpose) and isinstance(expr.operand, matchpy.Wildcard) and expr.operand.variable_name == name:
        return expr.operand
    elif isinstance(expr, InverseTranspose) and isinstance(expr.operand, matchpy.Wildcard) and expr.operand.variable_name == name:
        return Inverse(expr.operand)
    elif isinstance(expr, matchpy.Wildcard):
        return expr
    else:
        return type(expr)(*[remove_transpose(op, name) for op in expr.operands])
    
class OutputOperand():
    """docstring for OutputOperand"""
    def __init__(self, operand, storage_format):
        self.operand = operand
        self.storage_format = storage_format

    def __repr__(self):
        return "".join(["OutputOperand(", self.operand.name, ", ", self.storage_format.name, ")"])

class KernelVariant():
    """docstring for KernelVariant"""
    def __init__(self, arg_vals):
        self.arg_vals = arg_vals

class ExpressionKV(KernelVariant):
    """Contains information about expression variants.

    "expression variant" is used for BLAS kernel variants where the expression
    that is computed depends on an argument. For example, the SYMM kernel
    computes either C <- AB+C or C <- BA+C, depending on whether the "side"
    argument is "L" or "R".

    Attributes:
        arg_name (str): The name of the corresponding argument in the signature.
            If this is None, then there are no expression variants for  this
            kernel.
        arg_vals (dict): Maps argument values (str) to the corresponding
            expressions (Expression).
    """
    def __init__(self, arg_name, arg_vals):
        super().__init__(arg_vals)
        self.arg_name = arg_name
        
class PropertyKV(KernelVariant):
    """Contains information about property variants.

    "property variant" is used for BLAS kernel variants where the properties of
    some operands depend on an argument. For TRSM for example, A is either upper
    or lower triangular, depending on whether the "uplo" argument is "U" or
    "L".

    Note:
        Properties of symbols that do not depend on arguments (e.g. for SYMM, A
        is always symmetric) are specified by setting that property on the
        symbol in the expression for ExpressionKV.

    Attributes:
        arg_name (str): The name of the corresponding argument in the signature.
        operand (Expression): The object specifies this operand's properties.
        wildcard_name (str): The name of the wildcard corresponding to operand.
    """
    def __init__(self, arg_name, arg_vals, operand):
        super().__init__(arg_vals)
        self.arg_name = arg_name
        self.operand = operand
        self.wildcard_name = to_wildcard_name(operand)

class DefaultValueKV(KernelVariant):
    """Contains information about default value variants.

    "default value variant" is used for BLAS kernel variants where one of the
    operands is replaced with a default value. Which default values make sense
    depends on the kernel. For GEMM for example (C <- alpha AB + beta C), it
    is reasonable to use 1 for alpha and 0 and 1 for beta. Then, this kernel can
    also be used if there are no scalars and/or no addition. Otherwise, the
    kernel could only be used if every wildcard matches one operand. In all
    cases, one ReductionKernel is generated where no default values are used.

    Note:
        For each operand that is supposed to be replaced with default values,
        one separate DevaultValueKV is needed.

        If there are multiple DefaultValueKV objects, all possible combinations
        are generated.

    Attributes:
        arg_vals (dict): Maps argument values (Expression) to a list of
            corresponding default values.
        operand (Expression): The operand to be replaced with a default value.
        wildcard_name (str): The name of the wildcard corresponding to operand.
    """
    def __init__(self, operand, values):
        # The following entry guarantees that the case where no default value
        # is used is still considered.
        self.arg_vals = {None: None}
        for value in values:
            self.arg_vals[str(value.value)] = value
        self.operand = operand
        self.wildcard_name = to_wildcard_name(operand)

class OperatorKV(KernelVariant):
    """Contains information about operator variants.

    "operator variant" is used for BLAS kernels where there are multiple options
    for a unary operator. For GEMM for example, in C <- Op(A)*Op(A), Op can
    either be the identity function of transposition, depending on the values
    of transA and transB.

    Note:
        Instead of Op, Op1 and Op2 have to be used in the expression to identify
        the different operators.

    Attributes:
        arg_name (str): The name of the corresponding argument in the signature.
        operator (Operator): The name of the operator placeholder (either Op1 or
        Op2) to be replaced with a real operator.

    """
    def __init__(self, arg_name, arg_vals, operator):
        super().__init__(arg_vals)
        self.arg_name = arg_name
        self.operator = operator


############################
# New Expression classes

class Op1(Operator):
    """docstring for Op1"""

    name = "Op1"
    arity = matchpy.Arity.unary
    associative = False
    commutative = False
    one_identity = False
    infix = False

    def __init__(self,  *operands, variable_name=None):
        super().__init__(*operands, variable_name=variable_name)


class Op2(Operator):
    """docstring for Op2"""

    name = "Op2"
    arity = matchpy.Arity.unary
    associative = False
    commutative = False
    one_identity = False
    infix = False

    def __init__(self,  *operands, variable_name=None):
        super().__init__(*operands, variable_name=variable_name)


############################
# Kernel description classes

class KernelDescriptionError(Exception):
    pass        

class KernelDescription():
    """Description of BLAS-like kernels.

    This class contains a complete description of a BLAS-like kernel. It is
    used to generate ReductionKernel objects for specific instances of this
    kernel. In this context, an "instance" is a fixed combination of arguments,
    for example whether a certain operand is transposed or not, or whether a
    matrix is upper or lower triangular. The motivation is that each such
    instance has to be represented by a separate pattern, but the task of
    creating all those patterns should not be left to the user.

    Note:
        It is important that the names and identifiers used in the signature,
        the KernelVariant objects and the Arguments objects are consistent.
        This is not checked.

        Properties of symbols that do not depend on arguments (e.g. for SYMM, A
        is always symmetric) are specified by setting that property on the
        symbol in the expression for ExpressionKV.

    Attributes:
        signature (CodeTemplate): A string containing the signature of the
            kernel. All argument names have to be preceded by "$".
        expr_variant (ExpressionKV): The expression variants of this kernel.
        variants (KernelVariant): All other variants of this kernel, excluding
            the ExpressionKV object.
        arguments (Argument): The arguments of the kernel.
        return_value (OutputOperand): The operand that contains the result of the
            kernel.
        kernel_types (KernelType): The kernel types to be generated for this
            kernel.
        cost_function (function): A cost function for this kernel.
        wildcards (dict): Mapping operand names to WildcardSymbol objects. The
            wildcard names the the names of the operands preceded by an
            underscore "_".

    """
    # TODO rename kernel variants: KV as prefix
    def __init__(self, expr_variant, variants, input_operands, return_value, cost_function, pre_code, signature, post_code, arguments, kernel_types=[KernelType.identity]):
        self.expr_variant = expr_variant
        self.variants = variants
        self.input_operands = input_operands
        self.return_value = return_value
        self.kernel_types = kernel_types
        self.cost_function = cost_function

        self.pre_code = CodeTemplate(pre_code)
        self.signature = CodeTemplate(signature)
        self.post_code = CodeTemplate(post_code)

        self.arguments = arguments

        self.wildcards = dict()
        for expr in self.expr_variant.arg_vals.values():
            for node, _ in expr.preorder_iter():
                if isinstance(node, Symbol) and node.name not in self.wildcards:
                    # self.wildcards[node.name] = matchpy.Wildcard.symbol(to_wildcard_name(node))
                    self.wildcards[node.name] = matchpy.Wildcard.symbol(to_wildcard_name(node), symbol_type=type(node))
        
        self.wildcards[self.return_value.operand.name] = matchpy.Wildcard.symbol(to_wildcard_name(self.return_value.operand), symbol_type=type(self.return_value.operand))

        # Changing identifiers in the signature from symbol name to wildcard name.
        self.signature.substitute_identifiers({symbol_name: wildcard.variable_name for symbol_name, wildcard in self.wildcards.items()})

        # print(self.wildcards)
        # Replacing symbols in the Argument objects with wildcards.
        for arg in self.arguments:
            arg.operand = substitute_symbols_with_wildcards(arg.operand, self.wildcards)
            # arg.operand = self.wildcards[arg.operand.name]

        for input_operand in self.input_operands:
            input_operand.operand = self.wildcards[input_operand.operand.name]

        self.return_value.operand = self.wildcards[self.return_value.operand.name]

    def generate_kernels(self):
        """Generator for ReductionKernel objects.

        Yields:
            ReductionKernel: Each object represents one instance of a BLAS-like
                kernel.
        """
        # TODO using dictionaries for arg_vals does not make a lot of sense
        # the lookup functionality is never used!

        types = [type(variant) for variant in self.variants]
        # iterating over expression variants
        for arg, expr in self.expr_variant.arg_vals.items():

            arguments_copy1 = [copy.copy(arg) for arg in self.arguments] # copy is not enough, deepcopy is not necessary

            if self.expr_variant.arg_name is not None:
                # This kernel does have expression variants.
                # TODO I don't like the previous comment.
                arguments_copy1.append(ConstantArgument(self.expr_variant.arg_name, arg))

            # Dealing with the remaining kernel variants
            # Step 1: Generating all combinations.
            for arg_val_pairs in itertools.product(*[variant.arg_vals.items() for variant in self.variants]):
                expr_copy2 = expr
                arguments_copy2 = [copy.copy(arg) for arg in arguments_copy1] # copy is not enough, deepcopy is not necessary

                constraints_dict = dict()

                kernel_io = KernelIO()

                # Step 2: For each combination, iterating over all kernel variant objects.
                for i, _type in enumerate(types):
                    if _type == PropertyKV:
                        # Adding property constraints to constraint object
                        arg, property = arg_val_pairs[i]
                        # Properties can not be stored on the symbol, because
                        # the same Symbol instance is used for multiple
                        # generated kernels (the alternative would be to
                        # traverse expr_copy2, which contains a copy of the
                        # symbol).
                        constraints_dict.setdefault(self.variants[i].wildcard_name, set()).add(property)
                        arguments_copy2.append(ConstantArgument(self.variants[i].arg_name, arg))
                    elif _type == DefaultValueKV:
                        # replacing symbol with default value
                        arg, value = arg_val_pairs[i]
                        if arg is not None:
                            expr_copy2 = replace_symbol(expr_copy2, self.variants[i].operand, value)
                            kernel_io.add_input(self.variants[i].operand, value, storage_format.StorageFormat.full)
                    elif _type == OperatorKV:
                        # replacing operator placeholder with actual operator
                        arg, new_operator = arg_val_pairs[i]
                        expr_copy2 = replace_operator(expr_copy2, self.variants[i].operator, new_operator)
                        for argument in arguments_copy2:
                            if argument.operand:
                                argument.operand = replace_operator(argument.operand, self.variants[i].operator, new_operator)
                        arguments_copy2.append(ConstantArgument(self.variants[i].arg_name, arg))
                
                expr_copy2 = simplify(expr_copy2)

                remaining_operands = []
                seen_before = set()
                for node, _ in expr_copy2.preorder_iter():
                    if isinstance(node, Symbol) and node.name not in seen_before:
                        remaining_operands.append(node)

                # print(expr_copy2)
                # print(self.signature)
                # print(remaining_operands)
                # print(kernel_io)
                # print(self.return_value.operand, self.return_value.storage_format)
                # print(self.wildcards)

                # Add property constraints (including Matrix, Vector and Scalar)
                # This has to be done after simplifying to make sure that only
                # those operands are included that actually show up in the
                # expression.
                constraints_list = []
                for operand in remaining_operands:
                    _wildcard_name = self.wildcards[operand.name].variable_name
                    if operand.properties or _wildcard_name in constraints_dict:
                        property_set = operand.properties.union(constraints_dict.get(_wildcard_name, set()))
                        # Those properties can be removed because variables use symbol_type
                        property_set.difference_update((properties.MATRIX, properties.VECTOR, properties.SCALAR))
                        if property_set:
                            constraints_list.append(PropertyConstraint(_wildcard_name, property_set))

                remaining_wildcards = [self.wildcards[operand.name] for operand in remaining_operands]
                
                remaining_input_operands = []
                for input_operand in self.input_operands:
                    if input_operand.operand in remaining_wildcards:
                        remaining_input_operands.append(input_operand)

                # print([io.operand for io in remaining_input_operands])

                # Replace symbols with wildcards in the expression
                expr_copy2 = replace_symbols(expr_copy2, self.wildcards)

                # print(constraints_list)

                # print(self.wildcards)
                # print(repr(expr_copy2))
                # print(constraints)

                kernel_io.replace_variables(self.wildcards)

                for kernel_type in self.kernel_types:
                    yield ReductionKernel(matchpy.Pattern(expr_copy2, *constraints_list), remaining_input_operands, self.return_value, self.cost_function, self.pre_code, self.signature, self.post_code, arguments_copy2, kernel_type, kernel_io)        

############################
# Auxiliary functions

def replace_symbol(expr, old, new):
    if isinstance(expr, Symbol): # only symbols have names
        if expr.name == old.name: # TODO what about expr == old?
            return new
        else:
            return expr
    elif isinstance(expr, matchpy.Wildcard):
        return expr
    else:
        return type(expr)(*[replace_symbol(operand, old, new) for operand in expr.operands])

def replace_symbols(expr, mapping):
    if isinstance(expr, Symbol): # only symbols have names
        if expr.name in mapping: # TODO what about expr == old?
            return mapping[expr.name]
        else:
            return expr
    elif isinstance(expr, matchpy.Wildcard):
        return expr
    else:
        return type(expr)(*[replace_symbols(operand, mapping) for operand in expr.operands])

def replace_operator(expr, old, new):
    # print(expr)
    # print(type(expr), old, isinstance(expr, old))
    if isinstance(expr, old):
    # if type(expr) == old:
        return replace_operator(new(*expr.operands), old, new)
    elif not isinstance(expr, Symbol) and not isinstance(expr, matchpy.Wildcard):
        return type(expr)(*[replace_operator(operand, old, new) for operand in expr.operands])
    else:
        return expr

if __name__ == "__main__":
    
    from patternmatcher.functional import Pattern

    from patternmatcher.expression import Matrix, Transpose, Scalar, NumericConstant, simplify, Inverse

    import patternmatcher.InferenceOfProperties


    A = KMatrix("A")
    A.set_property(properties.SQUARE)
    B = KMatrix("B")
    alpha = KScalar("alpha")
    # one = NumericConstant(0)
    patt = Pattern(Times(A, B)) 
    cf = lambda md: 1

    # kernel = Kernel(patt, cf, B, KernelType.transpose)
    # print(kernel.operation_template)
    # print(kernel.pattern_with_context)
    # print(kernel.pattern)
    # print(kernel.replacement_with_context_template)
    # print(kernel.replacement_template)
    # print(kernel.overwriting)

    # print(replace_operator(Times(Op1(A), B), Op1, Transpose))
    # print(replace_symbol(Times(alpha, A, B), alpha, one))

    # gemm = KernelDescription("gemm($transA, $transB, $M, $N, $K, $alpha, $A, $ldA, $B, $ldB, $beta, $C, $ldC)"
    #                          [ExpressionKV(Times(A, B))],
    #                          [])

    trsm = KernelDescription("trsm($side, $uplo, $transA, $diag, $M, $N, $alpha, $A, $ldA, $B, $ldB)",
                             ExpressionKV(
                                "side",
                                {"L": Times(alpha, Op1(Inverse(A)), B),
                                 "R": Times(alpha, B, Op1(Inverse(A)))}
                                ),
                             [OperatorKV(
                                "transA",
                                {"N": Transpose,
                                 "T": Transpose},
                                Op1
                                ),
                              DefaultValueKV(
                                alpha,
                                [NumericConstant(1)]
                                ),
                              PropertyKV(
                                "uplo",
                                {"U": properties.UPPER_TRIANGULAR,
                                 "L": properties.LOWER_TRIANGULAR},
                                A
                                )
                             ],
                             [], # Argument objects
                             B, # return value
                             cf # cost function
                            )

    for var in trsm.generate_kernels():
        pass