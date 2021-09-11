from ...algebra.expression import Symbol

from ... import config

from ...code_generation.utils import MatchedKernel
from ...code_generation.memory import memory as memory_module
from ...code_generation.memory import storage_format as sf

import matchpy

class AbstractKernel():
    """docstring for AbstractKernel"""

    _counter = 0

    def __init__(self, cost_function, signature, arguments):
        self.id = AbstractKernel._counter
        AbstractKernel._counter += 1
        self._cost_function = cost_function

        self.signature = signature

        # list of Argument objects
        self.arguments = arguments

    def set_match(self, match_dict):

        matched_kernel = MatchedKernel()

        #############
        # Signature

        matched_kernel.signature = self.signature

        #############
        # Arguments

        matched_kernel.arguments = self.arguments

        #############
        # Cost

        matched_kernel.cost = self.cost(match_dict)

        return matched_kernel

    def cost(self, match_dict):
        """Returns the cost of this operation.
        """
        arg_dict = dict()
        for arg in self.arguments:
            if isinstance(arg, SizeArgument):
                arg_dict[arg.name] = arg.get_value(match_dict)
        return self._cost_function(arg_dict)



class PropertyConstraints():
    """docstring for PropertyConstraints
    
    Acts as a constraint 'function' for wildcards in kernels.
    """
    def __init__(self):
        self.property_dict = dict()

    def add_constraints(self, wildcard_name, *properties):
        self.property_dict.setdefault(wildcard_name, []).extend(properties)

    def __call__(self, match_dict):
        return all(all(match_dict[wildcard_name].has_property(prop) for prop in properties) for wildcard_name, properties in self.property_dict.items())

    def __repr__(self):
        return "PropertyConstraints(" + repr(self.property_dict) + ")" 

class InputOperand():
    """docstring for InputOperand"""
    def __init__(self, operand, storage_format):
        self.operand = operand
        self.storage_format = storage_format

    def __repr__(self):
        return "".join(["InputOperand(", self.operand.name, ", ", self.storage_format.name, ")"])

class Argument():
    """docstring for Argument"""
    def __init__(self, name, operand):
        self.name = name
        self.operand = operand

    def get_value(self, match_dict, memory):
        # TODO it looks like SizeArgument are the only ones that actually need
        # this funtion.
        raise NotImplementedError()

    def get_replacement(self, match_dict, memory, storage_formats):
        """

        Returns:
            str: A line of code that has to be placed before the kernel call.
            str: A replacement for the placeholder for this argument in the
                signature.
        """
        raise NotImplementedError()

class SizeArgument(Argument):
    """docstring for SizeArgument"""
    def __init__(self, name, operand, dimension, as_value=False):
        super().__init__(name, operand)
        self.dimension = dimension
        self.as_value = as_value

    def get_value(self, match_dict, memory=None):
        # size = self.operand.replace_copy(match_dict).size
        # TODO don't use substitute here
        op = matchpy.substitute(self.operand, match_dict)
        if self.dimension == "rows":
            return op.rows
        elif self.dimension == "columns":
            return op.columns
        elif self.dimension == "entries":
            return op.rows*op.columns
        else:
            raise ValueError("{} is not a valid dimension.".format(self.dimension))
        # return operator.attrgetter(self.dimension)(match_dict[self.operand.name])

    def get_replacement(self, match_dict, memory=None, storage_formats=None):
        if config.c:
            value = self.get_value(match_dict)
            if self.as_value:
                return None, str(value)
            else:            
                return to_c_variable_definition(value), "".join(["&", to_c_variable(value)])
        elif config.julia:
            value = self.get_value(match_dict)
            return None, str(value)
        else:
            raise config.LanguageOptionNotImplemented()

    def __repr__(self):
        return "".join(["SizeArgument(", self.name, ", ", repr(self.operand), ", ", repr(self.dimension), ")"])

class StrideArgument(Argument):
    """docstring for StrideArgument"""
    def __init__(self, name, operand, dimension, as_value=False):
        super().__init__(name, operand)
        self.dimension = dimension
        self.as_value = as_value

    def get_value(self, match_dict, memory):
        matched_operand = match_dict[self.operand.name]
        try:
            mem_loc = memory.lookup[matched_operand.name]
        except KeyError:
            raise memory_module.OperandNotInMemory()

        if self.dimension == "rows":
            return mem_loc.stride[0]
        elif self.dimension == "columns":
            return mem_loc.stride[1]
        else:
            raise ValueError("{} is not a valid dimension.".format(self.dimension))

    def get_replacement(self, match_dict, memory, storage_formats=None):
        if config.c:
            value = self.get_value(match_dict, memory)
            if self.as_value:
                return None, str(value)
            else:            
                return to_c_variable_definition(value), "".join(["&", to_c_variable(value)])
        if config.julia:
            value = self.get_value(match_dict, memory)
            return None, str(value)
        else:
            raise config.LanguageOptionNotImplemented()

    def __repr__(self):
        return "".join(["StrideArgument(", self.name, ", ", repr(self.operand), ", ", repr(self.dimension), ")"])

class PropertyArgument(Argument):
    """docstring for PropertyArgument"""
    def __init__(self, name, operand, property, values):
        super().__init__(name, operand)
        self.property = property
        self.values = values

    def get_value(self, match_dict, memory=None):
        return self.property

    def get_replacement(self, match_dict, memory=None, storage_formats=None):
        has_property = match_dict[self.operand.variable_name].has_property(self.property)
        replacement = None
        if has_property:
            replacement = self.values[0]
        else:
            replacement = self.values[1]
        if config.c:
            return None, "".join(["\"", replacement, "\""])
        elif config.julia:
            return None, "".join(["\'", replacement, "\'"])
        else:
            raise config.LanguageOptionNotImplemented()

    def __repr__(self):
        return "".join(["PropertyArgument(", self.name, ", ", repr(self.operand), ", ", repr(self.property), ", ", repr(self.values), ")"])

class StorageFormatArgument(Argument):
    def __init__(self, name, operand, storage_formats):
        super().__init__(name, operand)
        self.storage_formats = storage_formats

    def get_replacement(self, match_dict, memory, storage_formats):
        replacement = None
        matched_operand_name = match_dict[self.operand.variable_name].name

        storage_format = storage_formats[matched_operand_name]
        try:
            replacement = self.storage_formats[storage_format]
        except KeyError:
            raise sf.IncompatibleStorageFormats("Expected one of {}, got {}.".format([format.name for format in self.storage_formats.keys()] , storage_format))

        if config.c:
            return None, "".join(["\"", replacement, "\""])
        elif config.julia:
            return None, "".join(["\'", replacement, "\'"])
        else:
            raise config.LanguageOptionNotImplemented()

    def __repr__(self):
        return "".join(["StorageFormatArgument(", self.name, ", ", repr(self.operand), ", ", repr(self.storage_formats), ")"])

class ConstantArgument(Argument):
    """docstring for ConstantArgument"""
    def __init__(self, name, value):
        super().__init__(name, None)
        self.value = value

    def get_value(self, match_dict=None, memory=None):
        return self.value

    def get_replacement(self, match_dict=None, memory=None):
        if config.c:
            value = self.get_value()
            if isinstance(value, str):
                return None, "".join(["\"", value, "\""])
            else:
                return to_c_variable_definition(value), "".join(["&", to_c_variable(value)])
        if config.julia:
            value = self.get_value(match_dict, memory)
            if isinstance(value, str):
                return None, "".join(["\'", value, "\'"])
            else:
                # For most BLAS wrappers in Julia, it is important that
                # constants are not integers, but Float64. Thus, they have to
                # be printed as "1.0", not "1". This happens automatically if
                # they are of type float in Python, which has to be made sure
                # then creating ConstantArgument objects. Instead of always
                # converting to float here, we assume that they already
                # have the right type (to account for functions that actually
                # require integers). For this to work, it is probably necessary
                # to also make sure that ConstantScalars always have float
                # values. Should any of this ever end up being a problem, it
                # might be a good solution to introduce different
                # ConstantArgument objects or to introduce a type attribute.
                return None, str(value)
        else:
            raise config.LanguageOptionNotImplemented()

    def __repr__(self):
        return "".join(["ConstantArgument(", self.name, ", ", repr(self.operand), ", ", repr(self.value), ")"])

# TODO BandwidthArgument?
        

def to_wildcard_name(operand):
    return "_{}".format(operand.name)

def substitute_symbols_with_wildcards(expr, _substitution_dict=dict()):
    if isinstance(expr, Symbol):
        return _substitution_dict[expr.name]
    expr.operands = [substitute_symbols_with_wildcards(operand, _substitution_dict) for operand in expr.operands]
    return expr

def to_c_variable(value):
    variable_name = None
    if isinstance(value, float):
        value_string = str(value).replace(".", "p")
        variable_name = "".join([config.data_type_string, value_string])
    elif isinstance(value, int):
        variable_name = "".join(["int", str(value)])
    variable_name = variable_name.replace("-", "m")
    return variable_name

def to_c_variable_definition(value):
    type_name = None
    if isinstance(value, float):
        type_name = config.data_type_string
    elif isinstance(value, int):
        type_name = "int"
    return "{0} {1} = {2};\n".format(type_name, to_c_variable(value), value)


def le_property_sets(s1, s2):
    """ <= relation on sets of properties.

    s1 <= s2 if s1 is a more restricted set of properties than s2.     
    
    For every element in the "larger" set, there is a smaller or equal element
    in the "smaller" set.

    For every property in the more general property set, there is a more
    specific (or equal) property in the smaller set. The empty set is the most
    general set.

    Examples:
        set([Property.SQUARE, Property.DIAGONAL]) < set([Property.LOWER_TRIANGULAR])

    Note:
        This function implements a partial ordering. Property sets may not be
        comparable. As an example, consider set(Property.LOWER_TRIANGULAR) and
        set(Property.UPPER_TRIANGULAR).
    """
    return all(any(e1 < e2 or e1 == e2 for e1 in s1) for e2 in s2)


class PropertyTuple(tuple):
    """Tuple of sets of properties with partial ordering.
    
    This subclass of tuple should only be used for sets of properties. It offers
    comparision functions which implement an extension of the ordering defined
    by le_property_sets(s1, s2) to tuples of sets of properties. This order is
    also called "product order" (https://en.wikipedia.org/wiki/Product_order).

    Note:
        An example for the usecase of this kernel:
        Take A*B as an example, where A is lower triangular. This can be
        computed with both gemm and trmm. Pattern matching will find both
        kernels. PropertyTuple is used to represent the property constraints of
        those kernels. To quickly identify that gemm is unnecessarily general,
        the partial ordering on properties is extended to PropertyTuple and used
        to identify that trmm is strictly more specific than gemm.

    """
    def __new__(cls, elements):
        return super().__new__(cls, tuple(elements))

    # def __init__(self, kernel, renaming, elements):
    #     # super().__init__()
    #     self.kernel = kernel
    #     self.renaming = renaming

    def _content_str(self):
        set_strs = ["{" + ", ".join([prop.name for prop in p_set]) + "}" for p_set in self]
        return ", ".join(set_strs)        

    def __str__(self):
        return "PT(" + self._content_str() + ")"

    def __repr__(self):
        return "PropertyTuple({0})".format(self._content_str())
 
    def __lt__(self, other):
        """Extension of le_property_sets to tuples of property sets.

        Note:
            This implements a partial ordering.
        """
        if self == other:
            return False
        else:
            return all(le_property_sets(e1, e2) for e1, e2 in zip(self, other))

    def __le__(self, other):
        """Extension of le_property_sets to tuples of property sets.

        Note:
            This implements a partial ordering.
        """
        return all(le_property_sets(e1, e2) for e1, e2 in zip(self, other))

    def __ge__(self, other):
        """Extension of le_property_sets to tuples of property sets.

        Note:
            This implements a partial ordering.
        """
        return all(le_property_sets(e1, e2) for e1, e2 in zip(other, self))   
