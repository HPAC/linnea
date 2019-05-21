from .expression import Symbol, Matrix, Scalar, \
                        Equal, Plus, Times, Transpose, Inverse, \
                        BlockedExpression, Operator

from .utils import scalar_subexpressions

import matchpy

class SizeMismatchInProduct(Exception):
    pass

class SizeMismatchInSum(Exception):
    pass

class SizeMismatchInEquation(Exception):
    pass

def check_validity(expr):
    """Checks if an expression is valid

    It is checked whether the dimensions of sums and products match. If not, an
    exception is raised.

    [TODO] at the moment, this function does not consider conflicting properties
    (that should probably happen in a seperate function)

    [TODO] at the moment, blcoked expressions are not considered
    """

    if isinstance(expr, Symbol):
        return True

    if isinstance(expr, Equal):
        expr1, expr2 = expr.operands
        if expr1.size != expr2.size:
            raise SizeMismatchInEquation("%s = %s with sizes %s and %s" \
                                    % (expr1, expr2, expr1.size, \
                                       expr2.size))

        return check_validity(expr1) and check_validity(expr2)

    if isinstance(expr, Plus):
        prev_term = expr.operands[0] # previous term is only stored for printing
        size = prev_term.size
        
        # compare size of first term with consecutive terms
        for term in expr.operands[1:]:
            if size != term.size:
                raise SizeMismatchInSum("%s + %s with sizes %s and %s" \
                                        % (prev_term, term, prev_term.size, \
                                           term.size))
            prev_term = term

        return all(check_validity(term) for term in expr.operands)

    if isinstance(expr, Times):

        scalar_expressions = dict(scalar_subexpressions(expr))
        compare_with = 0 # for printing exception
        factors = expr.operands

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
                elif factor.size[1]==factors[compare_with].size[0]:
                    valid_dims = True
            # if current position is end of scalar expression, don't check if
            # dimension match
            elif n in scalar_expressions.values():
                valid_dims = True
            # normal check
            else:
                compare_with = n+1
                if factor.size[1]==factors[compare_with].size[0]:
                    valid_dims = True
            if not valid_dims:
                raise SizeMismatchInProduct("%s * %s with sizes %s and %s" \
                                            % (factor, factors[compare_with], \
                                            factor.size, \
                                            factors[compare_with].size))

        return all(check_validity(factor) for factor in factors)
    if isinstance(expr, Operator) and expr.arity is matchpy.Arity.unary:
        return check_validity(expr.operand)
    return False
