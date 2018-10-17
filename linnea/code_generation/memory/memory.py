
# from clak.equations import Equations

from .storage_format import StorageFormat, select_storage_format_conversion
from .storage_format_conversions import in_place_conversion_dict, out_of_place_conversion_dict

from ...algebra.expression import Symbol, Constant, ConstantScalar, IdentityMatrix, ZeroMatrix
from ...algebra.properties import Property as properties

import copy
import itertools
import operator

from ... import config

class OperandNotInMemory(Exception):
    pass

class MemoryOperation():
    pass

class FreeMemory(MemoryOperation):
    """docstring for FreeMemory"""
    def __init__(self, mem_loc):
        self.mem_loc = mem_loc

    def __repr__(self):
        return " ".join(["free", self.mem_loc.name])

    def code(self):
        if config.c:
            return "free({0});\n".format(self.mem_loc.name)
        else:
            raise config.LanguageOptionNotImplemented()


class AllocateMemory(MemoryOperation):
    """docstring for AllocateMemory"""
    def __init__(self, mem_loc):
        self.mem_loc = mem_loc

    def __repr__(self):
        return " ".join(["allocate", self.mem_loc.name, str(self.mem_loc.size)])

    def code(self):
        if config.c:
            size = "*".join(str(val) for val in self.mem_loc.size)
            return "{0}* {1} = ({0}*) malloc(sizeof({0})*{2});\n".format(config.data_type_string, self.mem_loc.name, size)
        elif config.julia:
            size = ", ".join(str(val) for val in self.mem_loc.size)
            return "{0} = Array{{{1}}}(undef, {2})\n".format(self.mem_loc.name, config.data_type_string, size)
        else:
            raise config.LanguageOptionNotImplemented()


class InitializeConstant(MemoryOperation):
    """docstring for AllocateMemory"""
    def __init__(self, mem_loc, storage_format):
        self.mem_loc = mem_loc
        self.storage_format = storage_format

    def __repr__(self):
        return " ".join(["initialize", self.mem_loc.name, str(self.mem_loc.content[0]), str(self.mem_loc.size)])

    def code(self):
        if config.c:
            raise config.LanguageOptionNotImplemented()
        elif config.julia:
            operand = self.mem_loc.content[0]
            initialization = ""

            if isinstance(operand, ConstantScalar):
                initialization = str(operand.value)
            elif isinstance(operand, IdentityMatrix):
                if self.storage_format <= StorageFormat.full:
                    initialization = "eye({0}, {1}, {2})".format(config.data_type_string, *operand.size)
                elif self.storage_format == StorageFormat.diagonal_vector:
                    initialization = "ones({0}, {1})".format(config.data_type_string, min(operand.size))
            elif isinstance(operand, ZeroMatrix):
                initialization = "zeros({0}, {1}, {2})".format(config.data_type_string, *operand.size)

            return "{0} = {1}\n".format(self.mem_loc.name, initialization)
        else:
            raise config.LanguageOptionNotImplemented()

class CopyOperand(MemoryOperation):
    """docstring for CopyOperand"""

    # TODO stride is not considered at all yet
    # TODO in the best case, the generated code takes matrix properties into
    # account (triangular, diagonal, symmetric)
    # TODO what if source and destination have different stride? This could
    # happen when working with partitioned operands. In principle,
    # I can choose the strid of the destination, because it is newly allocated.
    # However, it probably makes most sense to use stride 1. 
    def __init__(self, source, destination):
        self.destination = destination
        self.source = source

    def __repr__(self):
        return " ".join(["copy", str(self.source.content), "from", self.source.name, "to", self.destination.name])

    def code(self):
        size = "*".join(str(val) for val in self.destination.size)
        if config.c:
            return "memcpy(&{0}, &{1}, sizeof({2})*{3});\n".format(self.destination.name, self.source.name, config.data_type_string, size)
        elif config.julia:
            return "blascopy!({0}, {1}, 1, {2}, 1)\n".format(size, self.source.name, self.destination.name)
        else:
            raise config.LanguageOptionNotImplemented()


class ConvertStorageFormat(MemoryOperation):
    """docstring for ConvertStorageFormat"""
    def __init__(self, source, destination, conversion):
        self.destination = destination
        self.source = source
        self.conversion = conversion

    def __repr__(self):
        return " ".join(["convert", self.source.name, self.conversion.source_format.name, "to", self.destination.name, self.conversion.target_format.name])

    def code(self):
        if len(self.destination.size) == 1:
            # I think in this case, the size is currently not needed.
            # If it is needed, what name do we use? length? n?
            return self.conversion.code_template.safe_substitute_str(input=self.source.name, output=self.destination.name, type=config.data_type_string)
        elif len(self.destination.size) == 2:
            rows, columns = self.destination.size
            return self.conversion.code_template.safe_substitute_str(input=self.source.name, output=self.destination.name, type=config.data_type_string, m=rows, n=columns)
        


class ConvertStorageFormatInPlace(MemoryOperation):
    """docstring for ConvertStorageFormat"""
    def __init__(self, mem_loc, conversion):
        self.mem_loc = mem_loc
        self.conversion = conversion

    def __repr__(self):
        return " ".join(["convert", self.mem_loc.name, self.conversion.source_format.name, "to", self.conversion.target_format.name, "(in place)"])

    def code(self):
        rows, columns = self.mem_loc.size
        return self.conversion.code_template.safe_substitute_str(op=self.mem_loc.name, type=config.data_type_string, m=rows, n=columns)
        

class Memory():
    """docstring for Memory"""
    _counter = 0

    # TODO what about a dictionary for free memory locations (including partial locations)?
    # when generating code, I could still free locations as soon as they are not needed anymore
    # However, then I need to touch the code/intermediate represenation another time

    # TODO
    # For matrices that are stored as vectors (diagonal, permutation),
    # the actual operand size is not equal to the space occupied in memory.

    def __init__(self, eqns):
        self.id = Memory._counter
        Memory._counter += 1

        self.locations = []
        self.lookup = dict() # Operand name (str): MemoryLocation
        self.storage_format = dict() # Operand name (str): set of properties (StorageFormat)

        # Creating MemoryLocations for all operands in expression and
        # initializing self.lookup and self.references.

        for equation in eqns:
            for _expr, _ in equation.rhs.preorder_iter():
                if isinstance(_expr, Symbol) and not isinstance(_expr, Constant):
                    try:
                        mem_loc = self.lookup[_expr.name]
                    except KeyError:
                        mem_loc = MemoryLocation(_expr, StorageFormat.full)
                        mem_loc.content.append(_expr)
                        self.locations.append(mem_loc)
                        self.lookup[_expr.name] = mem_loc
                        # self.references[_expr.name] = 1
                        self.storage_format[_expr.name] = StorageFormat.full
                        # if _expr.has_property(properties.SYMMETRIC):
                        #     # TODO this is an assumption
                        #     self.storage_format[_expr.name].add(StorageFormat.lower_triangular)
                    else:
                        # self.references[_expr.name] += 1
                        pass
                        # TODO is it necessary to do something here?

    # TODO do I really need this?
    def has_property(self, operand, memory_property):
        return memory_property in self.properties[operand.name]

    def add_operation(self, kernel_io, liveness=None, line_number=None):
        """Updates the memory according to the input and output of an operation.

        This includes updating all references, creating new memory locations for
        new operands and deleting empty ones.

        Args:
            kernel_io (KernelIO)

        Returns:
            memory_operations_before (list):
                Contains MemoryOperation objects that represent operations that
                have to be executed before the added operation.
            memory_operations_after (list):
                Contains MemoryOperation objects that represent operations that
                have to be executed before the added operation.
            operand_mapping (dict):
                Maps argument names (str) to memory location names (str). To be
                used to replace argument names in the signature.
        """

        # print(input, output)

        memory_operations_before = []
        memory_operations_after = []
        # print(self.content_string_with_format())
        for operand, storage_format in kernel_io.input_operands:
            if not isinstance(operand, Constant):
                """
                Is it ok that storage conversion happens while the references are
                updated? At the moment, storage conversion does not touch the
                references. Should if ever to that, one should happen after the
                other.

                At the moment, this more or less works "by accident" if one kernel
                recieved the same operand multiple times as input in case both need
                the same storage format. For the second one,
                self.storage_format[operand.name] already is the required format, so
                only one conversion happens. Should one kernel ever take the same
                operands with different formats, this does not work (in that case,
                there is a somewhat serious problem because most of what is
                happening here is based on the assumption that each operand is in
                memory at most once).
                """
                if not storage_format <= self.storage_format[operand.name]:
                    # print("storage format conversion necessary")
                    # print(storage_format, self.storage_format[operand.name])
                    memory_operations_before.extend(self.storage_conversion(operand, storage_format))

        """
        Create mapping of operand names to memory location names. This has to be
        done after the storage format conversion operations because they change
        memory locations.

        This also has to be done in two parts. For the input operands before
        creating the memory operations, for the output operands after that.

        TODO: This should be called operand_argument_mapping or so, because it
        maps the names of those arguments that refer to operands (as opposed to
        arguments about operand sizes, stride, side, uplo,...) to memory
        locations.
        """
        operand_mapping = dict()
        for variable, (input, _output) in kernel_io.entries.items():
            if input:
                input_operand, storage_format = input
                if isinstance(input_operand, Constant):

                    if _output:
                        # If the constant is passed as an argument that is
                        # overwritten, it has to be placed in a new memory
                        # location.

                        new_location = MemoryLocation(input_operand, storage_format)
                        new_location.content.append(input_operand)
                        self.locations.append(new_location)

                        self.lookup[input_operand.name] = new_location

                        self.storage_format[input_operand.name] = copy.copy(storage_format)

                        # For C, we would need two memory operations here: Allocation and initialization.
                        memory_operations_before.append(InitializeConstant(copy.deepcopy(new_location), storage_format))

                        operand_mapping[variable.variable_name] = new_location.name

                    else:
                        # If the constant is not overwritten, it can be passed
                        # directly.
                        if config.julia:
                            if isinstance(input_operand, ConstantScalar):
                                operand_mapping[variable.variable_name] = str(input_operand.value)
                            elif isinstance(input_operand, IdentityMatrix):
                                if storage_format == StorageFormat.full:
                                    operand_mapping[variable.variable_name] = "eye({0}, {1}, {2})".format(config.data_type_string, *input_operand.size)
                                elif storage_format == StorageFormat.diagonal_vector:
                                    operand_mapping[variable.variable_name] = "eye({0}, {1})".format(config.data_type_string, min(input_operand.size))
                            elif isinstance(input_operand, ZeroMatrix):
                                operand_mapping[variable.variable_name] = "zeros({0}, {1}, {2})".format(config.data_type_string, *input_operand.size)
                        else:
                            raise LanguageOptionNotImplemented()
                else:
                    operand_mapping[variable.variable_name] = self.lookup[input_operand.name].name

        # print(operand_mapping)

        # group by overwriting
        _output = list(kernel_io.output_operands())
        _output.sort(key=keyfunc) # custom keyfunc is necessary to deal with None
        output_grouped = []
        for key, group in itertools.groupby(_output, operator.itemgetter(1)):
            output_grouped.append(([(operand, storage_format) for operand, _, storage_format in group], key))
        # print(output_grouped)

        for operands, overwriting in output_grouped:
            if overwriting is None:
                for operand, storage_format in operands:
                    # Nothing is overwritten.
                    # Create new MemoryLocation for each operand.

                    new_location = MemoryLocation(operand, storage_format)
                    new_location.content.append(operand)
                    self.locations.append(new_location)

                    self.lookup[operand.name] = new_location

                    if storage_format is StorageFormat.symmetric_triangular_out:
                        self.storage_format[operand.name] = StorageFormat.symmetric_lower_triangular
                    else:
                        self.storage_format[operand.name] = storage_format

                    # code generation
                    memory_operations_before.append(AllocateMemory(copy.deepcopy(new_location)))

            elif liveness[overwriting.name][1]==line_number:
                # Something is overwritten, but the operand is not used anymore.

                location = self.lookup[overwriting.name]
                location.content.remove(overwriting)
                
                for operand, storage_format in operands:

                    location.content.append(operand)
                    self.lookup[operand.name] = location
                    if storage_format is StorageFormat.as_overwritten:
                        self.storage_format[operand.name] = self.storage_format[overwriting.name]
                    elif storage_format is StorageFormat.symmetric_triangular_out and not self.storage_format[overwriting.name] < StorageFormat.symmetric_triangular_out:
                        self.storage_format[operand.name] = StorageFormat.symmetric_lower_triangular
                    else:
                        self.storage_format[operand.name] = storage_format

                # code generation: none

            else:
                # Something is overwritten, and the operand is still needed.
                # The overwritten operand is copied into a new MemoryLocation.
                # We overwrite the operand in the existing location.
                existing_location = self.lookup[overwriting.name]
                existing_location_copy = copy.deepcopy(existing_location)

                # Create new location.
                new_location = MemoryLocation.from_existing(existing_location)
                new_location_copy = copy.deepcopy(new_location)

                # Move overwritten operand to new location.
                new_location.content.append(overwriting)
                # This loop seems to be completely unnecessary.
                # for _operand in new_location.content:
                #     self.lookup[_operand.name] = new_location
                self.lookup[overwriting.name] = new_location
                self.locations.append(new_location)

                existing_location.content.remove(overwriting)
                
                # Add output (i.e. the new operands) to the existing location.
                for operand, storage_format in operands:
                    existing_location.content.append(operand)
                    self.lookup[operand.name] = existing_location

                    if storage_format is StorageFormat.as_overwritten:
                        self.storage_format[operand.name] = self.storage_format[overwriting.name]
                    elif storage_format is StorageFormat.symmetric_triangular_out and not self.storage_format[overwriting.name] < StorageFormat.symmetric_triangular_out:
                        self.storage_format[operand.name] = StorageFormat.symmetric_lower_triangular
                    else:
                        self.storage_format[operand.name] = storage_format

                # code generation
                memory_operations_before.append(AllocateMemory(new_location_copy))
                memory_operations_before.append(CopyOperand(existing_location_copy, new_location_copy))

        # print(input.operands)
        for operand, _ in kernel_io.input_operands:
            if not isinstance(operand, Constant) and liveness[operand.name][1]==line_number:
                mem_op = self._remove_operand(operand)
                if mem_op:
                    memory_operations_after.append(mem_op)

        """
        Adding output operands to operand_mapping. This has to be done after
        generation all the memory operations because it's possible that the
        memory locations of output operands don't exist before.
        """
        for variable, (input, _output) in kernel_io.entries.items():
            if not input and _output:
                # It is ok here to look at one (arbitrary) output operand. We
                # only care about the memory location here, and they all have to
                # be in the same one. (The reason is that this is what KernelIO
                # describes: Which operands are passed via which argument. This
                # raises the interesting question: Is it possible to construct
                # the operand_mapping in a different place, such that it's not 
                # necessary to go through self.lookup?)
                output_operand = _output[0][0]
                if not isinstance(output_operand, Constant):
                    operand_mapping[variable.variable_name] = self.lookup[output_operand.name].name
        # print(operand_mapping)

        return memory_operations_before, memory_operations_after, operand_mapping


    def _remove_operand(self, operand):
        """Revomes an operand from memory.

        If necessary, it also removes the memory location. If the operand is not
        in memory, nothing is done and None is returned.
        """
        memory_operation = None

        try:
            location = self.lookup[operand.name]
        except KeyError:
            return None
        
        # TODO is this test necessary?
        if operand in location.content:
            location.content.remove(operand)
        
        # TODO am I sure that I want to remove them? I could also reuse them
        if not location.content:
            # As long as the memory location is not reused, it's not necessary to copy it.
            if config.c:
                memory_operation = FreeMemory(location)
            self.locations.remove(location)

        del self.lookup[operand.name]
        # del self.references[operand.name]
        del self.storage_format[operand.name]

        return memory_operation

    def storage_conversion(self, operand, target_format):
        """Generates operations to convert the storage format of an operand.

        Generates all necessary memory operations to convert the storage format
        of operand. It is automatically determined whether this conversion can
        be done in place or not.

        Args:
            operand (Expression): The operand that needs to be converted.
            target_format (StorageFormat): The storage format that operand needs
                to be converted to.

        Returns:
            List of MemoryOperation objects.
        """

        memory_operations = []

        """
        TODO introduce different classes for in place and out of place conversion.
        The respective classes throw an exception if the required conversion is
        not available. Out of place is exclusively for factorization object (Am
        I sure? That's also the for diag and permutation)
        Select converion outside of of that class (i.e. here), also to figure out
        the target format. ConvertStorageFormat does nothing but looking up and
        generating code.

        Three cases
        - in place
        - out of place
        - copy
        """

        if config.julia and (
            self.storage_format[operand.name] < StorageFormat.factorization_obj
            or (
                self.storage_format[operand.name] < StorageFormat.as_vector or 
                target_format < StorageFormat.as_vector
                )
            ):
            # Those conversion can only be done out of place.

            existing_location = self.lookup[operand.name]
            existing_location_copy = copy.deepcopy(existing_location)

            # new_location = MemoryLocation.from_existing(existing_location)
            new_location = MemoryLocation(operand, target_format)
            new_location_copy = copy.deepcopy(new_location)

            new_location.content.append(operand)

            self.lookup[operand.name] = new_location
            self.locations.append(new_location)

            existing_location.content.remove(operand)

            # print(self.storage_format[operand.name], target_format)
            conversion = select_storage_format_conversion(self.storage_format[operand.name], target_format, out_of_place_conversion_dict)

            memory_operations.append(ConvertStorageFormat(existing_location_copy, new_location_copy, conversion))
            if config.c and not existing_location.content:
                memory_operations.append(FreeMemory(existing_location_copy))
        elif (len(self.lookup[operand.name].content) == 1):
            # No other operand is in this location, so the conversion can be
            # done in place.
            """
            TODO
            Improve this. Eventually, we should look at whether the conversion
            would conflict with other operands stored in the same location. Do
            in place conversion if it doesn't.
            """
            # print(operand, self.storage_format[operand.name], target_format)
            conversion = select_storage_format_conversion(self.storage_format[operand.name], target_format, in_place_conversion_dict)

            mem_loc = self.lookup[operand.name]
            mem_loc_copy = copy.deepcopy(mem_loc)
            memory_operations.append(ConvertStorageFormatInPlace(mem_loc_copy, conversion))
        else:
            # Operand is copied to a new location and the conversion is
            # performed on the new location.
            # This is done if there are multiple operands in the same location.
            existing_location = self.lookup[operand.name]
            existing_location_copy = copy.deepcopy(existing_location)

            new_location = MemoryLocation.from_existing(existing_location)
            new_location_copy = copy.deepcopy(new_location)

            new_location.content.append(operand)

            self.lookup[operand.name] = new_location
            self.locations.append(new_location)

            conversion = select_storage_format_conversion(self.storage_format[operand.name], target_format, in_place_conversion_dict)

            existing_location.content.remove(operand)
            memory_operations.append(AllocateMemory(new_location_copy))
            memory_operations.append(CopyOperand(existing_location_copy, new_location_copy))
            memory_operations.append(ConvertStorageFormatInPlace(new_location_copy, conversion))
            # It is never necessary to free here, since this type of conversion
            # only happens is there is something else in the same memory
            # location.

        self.storage_format[operand.name] = conversion.target_format

        return memory_operations
        

    def update_references(self, eqns):
        """Updates the references.

        This function updates the references with the reference counts for eqns.
        If the number of references for an operand reaches zero, it is removed.

        The purpose of this function is to take into account that some
        derivation steps, such as common subexpression elminiation, change the
        number of occurences of some operands, even though nothing is computed.
        """
        memory_operations = []

        new_references = dict()
        for equation in eqns:
            for _expr, _ in equation.rhs.preorder_iter():
                if isinstance(_expr, Symbol) and not isinstance(_expr, Constant):
                    try:
                        refs = new_references[_expr.name]
                    except KeyError:
                        new_references[_expr.name] = 1
                    else:
                        new_references[_expr.name] += 1

        # go through old dict
        # in new dict, check number of references
        # if it is zero, remove operand
        # if location is empty, remove location
        for operand_name in self.references:
            if operand_name not in new_references:
                mem_loc = self.lookup[operand_name]
                operand = None
                idx = None
                for _idx, _operand in enumerate(mem_loc.content):
                    if operand_name == _operand.name:
                        # TODO this could be simplified: There is no need to
                        # for the assignment below, one could directly use
                        # _operand outside the loop. Is this considered bad style?
                        idx = _idx
                        break

                del mem_loc.content[idx]

                if not mem_loc.content:
                    self.locations.remove(mem_loc)
                    # As long as the memory location is not reused, it's not necessary to copy it.
                    if config.c:
                        memory_operations.append(FreeMemory(mem_loc))

                del self.lookup[operand_name]

        self.references = new_references

        return memory_operations


    def __str__(self):
        _lookup =  ", ".join([operand + ": " + mem_loc.name for operand, mem_loc in self.lookup.items()])
        _str = "".join([_lookup, "\n"])
        _storage = ""
        for op, format in self.storage_format.items():
            _storage_entry = str(op) + ": " + str(format.name) + " "
            _storage = "".join([_storage, _storage_entry])
        _str = "".join([_str, _storage, "\n"])
        _locations = "\n".join([str(mem_loc) + " [" + ", ".join([str(self.references[operand.name]) for operand in mem_loc.content ]) + "]" for mem_loc in self.locations])
        _str = "".join([_str, _locations]) 
        return _str

    def content_string(self):
        return ", ".join(["{}: {}".format(operand, mem_loc.name) for operand, mem_loc in self.lookup.items()])
        
    def content_string_with_format(self):
        return ", ".join(["{}: {}, {}".format(operand, mem_loc.name, self.storage_format[operand].name) for operand, mem_loc in self.lookup.items()])

def keyfunc(_tuple):
    item1 = _tuple[1]
    if not item1:
        return ""
    return item1.name

class MemoryLocation():
    """docstring for MemoryLocation"""

    _counter = 0

    def __init__(self, operand, storage_format):
        self.id = MemoryLocation._counter
        MemoryLocation._counter += 1
        self.name = "".join(["ml", str(self.id)])
        if storage_format < StorageFormat.as_vector:
            self.size = (min(operand.size),) # this has to be a tuple
            self.stride = (1,) # same here
        elif operand.has_property(properties.VECTOR):
            self.size = (max(operand.size),) # this has to be a tuple
            # max because operand.size is something like (10, 1)
            self.stride = (1,) # same here
        else:
            size = operand.size
            self.size = size
            self.stride = (size[0], 1)
        self.content = []

    @staticmethod
    def from_existing(existing_location):
        new_location = object.__new__(MemoryLocation)
        new_location.id = MemoryLocation._counter
        MemoryLocation._counter += 1
        new_location.name = "".join(["ml", str(new_location.id)])
        new_location.size = existing_location.size
        new_location.stride = existing_location.stride
        new_location.content = []
        return new_location

    def __str__(self):
        return " ".join([self.name, str(self.size), str(self.content)])

    def __deepcopy__(self, memo):
        cpy = object.__new__(type(self))
        cpy.id = self.id
        cpy.name = self.name
        cpy.content = self.content.copy()
        cpy.size = self.size
        cpy.stride = self.stride
        return cpy

# """
# Reminder: The properties in MemoryInput and MemoryOutput have to come from the
# description of the kernels.
# """

# class MemoryInput():
#     """Input to a kernel as needed by Memory.add_operation().

#     Describes the input of a kernel, containing all information necessary to
#     handle memory menagement (as done by the Memory class).

#     Note:
#         Elements in the same positions of the attribute lists belong together,
#         for example, storage_formats[0] is the required format for operands[0].

#     Attributes:
#         operands (list): Contains the input operands (Expression) of the kernel.
#         storage_formats (list): Contains the storage formats (StorageFormat)
#             required for the corresponding input operands.
#     """
#     def __init__(self):
#         self.operands = []
#         self.storage_formats = []

#     def add_operand(self, operand, storage_formats=None):
#         self.operands.append(operand)
#         self.storage_formats.append(storage_formats)

#     def __repr__(self):
#         content_str = ", ".join([str(c) for c in zip(self.operands, self.storage_formats)])
#         return "".join(["MemoryInput(", content_str, ")"])


# class MemoryOutput():
#     """Output of a kernel as needed by Memory.add_operation().

#     Describes the output of a kernel, containing all information necessary to
#     handle memory menagement (as done by the Memory class).

#     Note:
#         Elements in the same positions of the attribute lists belong together,
#         for example, storage_formats[0] is the output format of operands[0].

#     Attributes:
#         operands (list): Contains the output operands (Expression) of the
#             kernel.
#         overwriting (list): If the corresponding output operand overwrites one
#             of the input operands of the kernel, this input operand is stored
#             in this list. If nothing is overwritten, the value is None.
#         storage_formats (list): Contains the storage formats (StorageFormat)
#             of the output operands.
#     """
#     def __init__(self):
#         self.operands = []
#         self.overwriting = [] # expression or None
#         self.storage_formats = []

#     def add_operand(self, operand, overwriting, storage_formats=None):
#         self.operands.append(operand)
#         self.overwriting.append(overwriting)
#         self.storage_formats.append(storage_formats)

#     def __repr__(self):
#         content_str = ", ".join([str(c) for c in zip(self.operands, self.overwriting, self.storage_formats)])
#         return "".join(["MemoryOuput(", content_str, ")"])

class StorageDescription():
    """docstring for StorageDescription"""
    def __init__(self, *description):
        self.operands = []
        self.storage_formats = []
        if description:
            self.operands, self.storage_formats = zip(*description)

if __name__ == "__main__":

    import patternmatcher.InferenceOfProperties
    import patternmatcher.properties as properties
    from patternmatcher.expression import Matrix, Times, Plus, Transpose, Equal

    n = 10
    A = Matrix("A", (n, n))
    B = Matrix("B", (n, n))
    # B.set_property(properties.SYMMETRIC)
    C = Matrix("C", (n, n))
    D = Matrix("D", (n, n))
    E = Matrix("E", (n, n))
    # E.set_property(properties.SYMMETRIC)
    tmp1 = Matrix("tmp1", (n, n))
    L = Matrix("L", (n, n))
    U = Matrix("U", (n, n))
    X = Matrix("X", (n, n))
    eqns = Equations(Equal(X, Times(A, B, C, D, A)))
    # eqns = Equations(Equal(X, Times(A, B, C, D)))
    eqns2 = Equations(Equal(X, Times(A, A, B, C, D)))

    mem = Memory(eqns)
    
    
    mem_ops = []
    # mem_ops = mem.update_references(eqns2)

    mem.storage_format["A"] = StorageFormat.lower_triangular
    mem.storage_format["B"] = StorageFormat.symmetric_lower_triangular
    mem.storage_format["C"] = StorageFormat.full
    mem.storage_format["D"] = StorageFormat.full

    print(mem)

    # E <- A^T A + B
    input0 = MemoryInput()
    input0.add_operand(A, StorageFormat.full)
    input0.add_operand(A, StorageFormat.full)
    input0.add_operand(B, StorageFormat.symmetric_triangular)
    output0 = MemoryOutput()
    output0.add_operand(E, B, StorageFormat.as_overwritten)
    mem_ops.extend(mem.add_operation(input0, output0))

    # [A, B], [(tmp1, A)]
    input1 = MemoryInput()
    input1.add_operand(A, StorageFormat.full)
    input1.add_operand(B, StorageFormat.full)
    output1 = MemoryOutput()
    output1.add_operand(tmp1, A, StorageFormat.full)
    # mem_ops.extend(mem.add_operation(input1, output1))

    # [A, B], [(tmp1, B)]
    input2 = MemoryInput()
    input2.add_operand(A)
    input2.add_operand(B)
    output2 = MemoryOutput()
    output2.add_operand(tmp1, B)
    # mem_ops.extend(mem.add_operation(input2, output2))

    # [A, B], [(tmp1, None)]
    input3 = MemoryInput()
    input3.add_operand(A)
    input3.add_operand(B)
    output3 = MemoryOutput()
    output3.add_operand(tmp1, None)
    # mem_ops.extend(mem.add_operation(input3, output3))
    
    # [A], [(L, A), (U, A)]
    input4 = MemoryInput()
    input4.add_operand(A)
    output4 = MemoryOutput()
    output4.add_operand(L, A)
    output4.add_operand(U, A)
    # mem_ops.extend(mem.add_operation(input4, output4))

    # [L], [(tmp1, L)]
    input5 = MemoryInput()
    input5.add_operand(L)
    output5 = MemoryOutput()
    output5.add_operand(tmp1, L)
    # mem_ops.extend(mem.add_operation(input5, output5))

    # [L], [(E, L)]
    input6 = MemoryInput()
    input6.add_operand(L)
    output6 = MemoryOutput()
    output6.add_operand(E, L)
    # mem_ops.extend(mem.add_operation(input6, output6))

    # [L], [(E, None)]
    input7 = MemoryInput()
    input7.add_operand(L)
    output7 = MemoryOutput()
    output7.add_operand(E, None)
    # mem_ops.extend(mem.add_operation(input7, output7))

    # mem_ops.extend(mem.add_operation([D], [(tmp1, L)]))

    print(mem_ops)
    print(mem)
