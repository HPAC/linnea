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
