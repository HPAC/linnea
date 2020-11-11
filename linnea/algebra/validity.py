from .expression import Symbol, Matrix, Vector, Scalar, Equal, Plus, Times, \
                        Transpose, Inverse, InverseTranspose, Operator, \
                        Expression

from .properties import Property, negative_implications

from .utils import scalar_subexpressions

from ..frontend.export import export_expression

import matchpy

class ExpressionException(Exception):
    def __init__(self, message, *args):
        self.message = message
        self.args = args

    def __str__(self):
        return self.message.format(*self.args)

    def replace_expressions(self):
        new_args = []
        for arg in self.args:
            if isinstance(arg, Expression):
                arg = export_expression(arg, dict())
            new_args.append(arg)
        return type(self)(self.message, *new_args)

class ConflictingProperties(Exception):
    pass

class SizeMismatch(ExpressionException):
    pass

class InvalidExpression(ExpressionException):
    pass

def check_validity(expr, dependent_dimensions=False):
    """Checks if an expression is valid.

    It is checked whether the dimensions in sums, products and equations match.
    If not, an exception is raised.

    Args:
        dependent_dimensions (bool, optinal): If True, disables some checks
            that are not applicable if the input is only used to compute
            dependent dimensions.

    Returns:
        bool: True if the expression is valid, raises an exception otherwise.

    Raises:
        SizeMismatch: If operand sizes in any expression do not match.
        ConflictingProperties: If operands have conflicting properties.
    """

    if isinstance(expr, Symbol):
        properties = expr.properties
        false_properties = expr.false_properties

        if isinstance(expr, Vector):
            rows, columns = expr.size
            if rows == 1 and columns == 1:
                msg = "Vector {} has length 1."
                raise InvalidExpression(msg, expr)
        elif isinstance(expr, Matrix):
            rows, columns = expr.size
            if rows == 1:
                msg = "Matrix {} has one row. Use Vector instead."
                raise InvalidExpression(msg, expr)
            if columns == 1:
                msg = "Matrix {} has one column. Use Vector instead."
                raise InvalidExpression(msg, expr)          

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
        lhs, rhs = expr.operands
        if lhs.size != rhs.size:
            msg = "Size mismatch in equation: {} = {} with sizes {} and {}."
            raise SizeMismatch(msg, lhs, rhs, lhs.size, rhs.size)
        if not isinstance(lhs, Symbol):
            msg = "The left-hand side of an assignment cannot be an expression: {}."
            raise InvalidExpression(msg, expr)

        return check_validity(lhs, dependent_dimensions) and check_validity(rhs, dependent_dimensions)

    if isinstance(expr, Plus):
        prev_term = expr.operands[0] # previous term is only stored for printing
        size = prev_term.size
        
        # TODO wouldn't it be simpler to do this with a sliding window?
        # compare size of first term with consecutive terms
        for term in expr.operands[1:]:
            if size != term.size:
                msg = "Size mismatch in sum: {} + {} with sizes {} and {}."
                raise SizeMismatch(msg, prev_term, term, prev_term.size, term.size)
            prev_term = term

        return all(check_validity(term, dependent_dimensions) for term in expr.operands)

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
                msg = "Size mismatch in product: {} * {} with sizes {} and {}."
                raise SizeMismatch(msg, factor, factors[compare_with], factor.size, factors[compare_with].size)

        return all(check_validity(factor, dependent_dimensions) for factor in factors)
    if isinstance(expr, (Inverse, InverseTranspose)):
        if not expr.has_property(Property.SQUARE):
            msg = "Only square expressions can be inverted: {}."
            raise InvalidExpression(msg, expr)
        if not dependent_dimensions and not (expr.has_property(Property.FULL_RANK) or expr.has_property(Property.SPSD) or expr.has_property(Property.SCALAR)):
            msg = "Singular expressions cannot be inverted: {}."
            raise InvalidExpression(msg, expr)
        return check_validity(expr.operand, dependent_dimensions)
    if isinstance(expr, Operator) and expr.arity is matchpy.Arity.unary:
        return check_validity(expr.operand, dependent_dimensions)
    return False
