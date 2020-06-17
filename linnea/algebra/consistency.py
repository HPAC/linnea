from .properties import Property, negative_implications

from .expression import Equal, Symbol


class ConflictingProperties(Exception):
    pass

class ConflictingIndices(Exception):
    pass

def check_consistency(expression):
    if isinstance(expression, Equal):
        if expression.rhs.indices != expression.lhs.indices:
            msg = "Conflicting indices in {}".format(expression)
            raise ConflictingIndices(msg)

    for expr, _ in expression.preorder_iter():
        if isinstance(expr, Symbol):
            """Properties SQUARE, ROW_PANEL, and COLUMN_PANEL are set here to
            make sure that certain types of errors (for example rectangular SPD
            matrices) are detected.
            The more elegant solution would be to set those properties in
            Matrix.__init__(). However, this is currently not possible because
            the properties of matrices are used in the kernel description to
            generate the constraints for variables. There are many kernels that
            do not require a specific shape, but since we are using actual
            Matrix objects for the kernel description, they do have exactly one
            of those three shapes, leading to unwanted constraints.
            To fix this, the way KernelDescription objects work needs to be
            changed.
            Checking if properties are consistent could further be improved by
            computing the so called "deductive closure" of property sets, that
            is, adding all properties that are implied by others. The necessary
            implications currently do not exist for sets of properties
            (for example "diagonal and SPD implies symmetric").
            """
            rows, columns = expr.size
            if rows < columns:
                expr.set_property(Property.ROW_PANEL)
            elif rows > columns:
                expr.set_property(Property.COLUMN_PANEL)
            elif rows != 1: # Scalars must not be square.
                expr.set_property(Property.SQUARE)
            properties = expr.properties
            false_properties = expr.false_properties

            for prop in properties:
                intersection = properties.intersection(negative_implications[prop])
                if intersection:
                    msg = "{} has conflicting properties: {} contradicts {}.".format(expr, prop.value, ", ".join(p.value for p in intersection))
                    raise ConflictingProperties(msg)

            intersection = properties.intersection(false_properties)
            if intersection:
                msg = "{} has conflicting properties {}.".format(expr, ", ".join(p.value for p in intersection))
                raise ConflictingProperties(msg)
    return True
