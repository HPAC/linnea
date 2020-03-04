from .properties import negative_implications

from .expression import Equal, Symbol


class ConflictingProperties(Exception):
    pass

class ConflictingIndices(Exception):
    pass

def check_consistency(expr):

    if isinstance(expr, Equal):
        if expr.rhs.indices != expr.lhs.indices:
            msg = "Conflicting indices in %s" % (repr(expr),)
            raise ConflictingIndices(msg)

    for node, _ in expr.preorder_iter():
        if isinstance(node, Symbol):
            # properties, false_properties = TOS[node.name]
            properties = node.properties
            false_properties = node.false_properties
            intersection = properties.intersection(false_properties)
            if intersection:
                msg = "Property sets of %s have non-empty intersection: %s" % (repr(node), intersection)
                raise ConflictingProperties(msg)
            # print(node, properties)
            for prop in properties:
                intersection = properties.intersection(negative_implications[prop])
                if intersection:
                    msg = "%s has conflicting properties: %s contradicts %s" % (repr(node), prop, intersection)
                    raise ConflictingProperties(msg)
