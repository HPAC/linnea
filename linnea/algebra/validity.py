from .expression import Symbol, Matrix, Scalar, Equal, Plus, Times, Transpose, \
                        Inverse, InverseTranspose, Operator

from .properties import Property, negative_implications

from .utils import scalar_subexpressions

import matchpy

class ConflictingProperties(Exception):
    pass


class SizeMismatch(Exception):
    def __init__(self, message, error_info=None):
        super().__init__(message)

        self.error_info = error_info

class InvalidExpression(Exception):
    def __init__(self, message, error_info=None):
        super().__init__(message)

        self.error_info = error_info


def check_validity(expr):
    """Checks if an expression is valid.

    It is checked whether the dimensions in sums, products and equations match.
    If not, an exception is raised.

    Returns:
        bool: True if the expression is valid, raises an exception otherwise.

    Raises:
        SizeMismatch: If operand sizes in any expression do not match.
        ConflictingProperties: If operands have conflicting properties.
    """

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

    if isinstance(expr, Equal):
        expr1, expr2 = expr.operands
        if expr1.size != expr2.size:
            msg = "Size mismatch in equation: {} = {} with sizes {} and {}.".format(expr1, expr2, expr1.size, expr2.size)
            raise SizeMismatch(msg, ("assignment", expr1, expr2))

        return check_validity(expr1) and check_validity(expr2)

    if isinstance(expr, Plus):
        prev_term = expr.operands[0] # previous term is only stored for printing
        size = prev_term.size
        
        # TODO wouldn't it be simpler to do this with a sliding window?
        # compare size of first term with consecutive terms
        for term in expr.operands[1:]:
            if size != term.size:
                msg = "Size mismatch in sum: {} + {} with sizes {} and {}.".format(prev_term, term, prev_term.size, term.size)
                raise SizeMismatch(msg, ("sum", prev_term, term))
            prev_term = term

        return all(check_validity(term) for term in expr.operands)

    if isinstance(expr, Times):

        scalar_expressions = dict(scalar_subexpressions(expr))
        compare_with = 0 # for printing exception
        factors = expr.operands

        # TODO wouldn't it be simpler to do this with a sliding window?
        # compare number of columns of current factor with rows of next factor
        for n, factor in enumerate(factors[:-1]):
            valid_dims = False
            #print(str(n) + " " + str(factor))
            # if next position is beginning of scalar, skip to end of scalar to
            # check dimension
            if n+1 in scalar_expressions:
                compare_with = scalar_expressions[n+1]+1
                # if scalar at the end of product, don't check
                if compare_with == len(factors):
                    valid_dims = True
                elif factor.columns == factors[compare_with].rows:
                    valid_dims = True
            # if current position is end of scalar expression, don't check if
            # dimension match
            elif n in scalar_expressions.values():
                valid_dims = True
            # normal check
            else:
                compare_with = n+1
                if factor.columns == factors[compare_with].rows:
                    valid_dims = True
            if not valid_dims:
                msg = "Size mismatch in product: {} * {} with sizes {} and {}.".format(factor, factors[compare_with], factor.size, factors[compare_with].size)
                raise SizeMismatch(msg, ("product", factor, factors[compare_with]))

        return all(check_validity(factor) for factor in factors)
    if isinstance(expr, (Inverse, InverseTranspose)):
        if not expr.has_property(Property.SQUARE):
            msg = "Only square expressions can be inverted: {}.".format(expr)
            raise InvalidExpression(msg, expr)
        return check_validity(expr.operand)
    if isinstance(expr, Operator) and expr.arity is matchpy.Arity.unary:
        return check_validity(expr.operand)
    return False
