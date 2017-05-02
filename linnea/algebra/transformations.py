from . import expression as ae
from .properties import Property as properties

import operator
import matchpy
import itertools

class TransformationError(Exception):
    pass        

# @profile
def simplify(expr):
    if not isinstance(expr, ae.Symbol) and not isinstance(expr, matchpy.Wildcard):
        expr = type(expr)(*list(map(simplify, expr.operands)))
    if isinstance(expr, ae.Operator) and expr.arity is matchpy.Arity.unary:
        if isinstance(expr, ae.Identity):
            return expr.operand
        elif isinstance(expr, ae.Transpose):
            return transpose(expr.operand)
        elif isinstance(expr, ae.Inverse):
            return invert(expr.operand)
        elif isinstance(expr, ae.InverseTranspose):
            return invert_transpose(expr.operand)
        elif isinstance(expr, ae.Conjugate):
            return conjugate(expr.operand)
        elif isinstance(expr, ae.ConjugateTranspose):
            return conjugate_transpose(expr.operand)
        elif isinstance(expr, ae.InverseConjugate):
            return invert_conjugate(expr.operand)
        elif isinstance(expr, ae.InverseConjugateTranspose):
            return invert_conjugate_transpose(expr.operand)
    elif isinstance(expr, ae.Times):

        ############################
        # Product is zero if there is one zero factor.

        for operand in expr.operands:
            if operand.has_property(properties.ZERO) or (isinstance(operand, ae.ConstantScalar) and operand.value == 0):
                return ae.ZeroMatrix(*expr.size)

        scalars, non_scalars = expr.split_operands()

        ############################
        # Remove identity matrix

        if len(non_scalars) > 1:
            non_scalars = [operand for operand in non_scalars if not operand.has_property(properties.IDENTITY)]

        ############################
        # Computing product of ConstantScalars, sorting symbolic scalars.

        # TODO once expr^expr is possible, this should also compute compute
        # products symbolic scalars (if possible)
        if len(scalars) > 1:
            numeric_scalars = [s for s in scalars if isinstance(s, ae.ConstantScalar)]
            symbolic_scalars = [s for s in scalars if not isinstance(s, ae.ConstantScalar)]
            symbolic_scalars.sort()

            if len(numeric_scalars) > 1:
                new_value = 1
                for i in range(len(numeric_scalars)):
                    new_value *= numeric_scalars[i].value
                # itertools.accumulate(numeric_scalars, lambda x, y: x.value * y.value)
                numeric_scalars = [ae.ConstantScalar(new_value)]

            scalars = numeric_scalars + symbolic_scalars

        # remove scalar 1
        # It is sufficient to look at the first entry of scalars because if
        # there was more than one numeric scalar, they would have been
        # combined.
        if scalars and isinstance(scalars[0], ae.ConstantScalar) and scalars[0].value == 1:
            del scalars[0]
        
        ############################
        # Removing Q^T Q where Q has OrthogonalColumns property
        #
        # TODO equivalent thing for orthogonal rows

        i = 0;
        while i < len(non_scalars):
            if non_scalars[i].has_property(properties.ORTHOGONAL_COLUMNS):
                if i > 0 and isinstance(non_scalars[i-1], ae.Transpose) and non_scalars[i-1].operand == non_scalars[i]:
                    del non_scalars[i]
                    del non_scalars[i-1]
                    i -= 1
                    if i >= len(non_scalars):
                        break
            i += 1

        if len(non_scalars) == 0 and not expr.has_property(properties.SCALAR):
            non_scalars.append(ae.IdentityMatrix(*expr.size))

        ############################
        # Removing expr * expr^{-1}
        #
        # To minimize the number of expressions that have to be inverted,
        # only every second factor is inverted. So in a product
        # Times([A, B, C, D]), only B and D are inverted, and then compared
        # to the previous and next factors.

        if len(non_scalars) > 1:
            while True:
                l = len(non_scalars)
                # l-1 is used because for l = 3, we don't need i == 2 (i.e.
                # the last position). In that case, range = [0], and for 
                # i == 0, A1 == A2^{-1} and A2^{-1} == A3 is checked.
                for i in range(0, l-1, 2):
                    if non_scalars[i].inverse_of(non_scalars[i+1]):
                        if l > 2:
                            del non_scalars[i:i+2]
                        else:
                            size = (non_scalars[i].size[0], non_scalars[i+1].size[1])
                            non_scalars[i:i+2] = [ae.IdentityMatrix(*size)]
                            pass
                        # since indices change after removing something, we have
                        # to start over
                        break
                    elif i+2 < l and non_scalars[i+1].inverse_of(non_scalars[i+2]):
                        if l > 2:
                            del non_scalars[i+1:i+3]
                        else:
                            size = (non_scalars[i+1].size[0], non_scalars[i+2].size[1])
                            non_scalars[i+1:i+3] = [ae.IdentityMatrix(*size)]
                            pass
                        # since indices change after removing something, we have
                        # to start over
                        break
                else:
                    break

        operands = scalars + non_scalars
        return ae.Times(*operands)
    elif isinstance(expr, ae.Plus):

        if expr.has_property(properties.SCALAR):

            # print("before", expr)

            ############################
            # Dealing with scalar sums.
            #
            # If possible, numeric and symbolic scalars are added up here.
            
            terms = []
            new_value = 0
            for operand in expr.operands:
                if isinstance(operand, ae.Times):
                    if isinstance(operand.operands[0], ae.ConstantScalar):
                        # It is safe to assume that only the first operand can
                        # be a ConstantScalar because because the operands
                        # of the current expression were already simplified.
                        terms.append([operand.operands[0], operand.operands[1:]])
                    else:
                        terms.append([1, operand.operands])
                elif isinstance(operand, ae.ConstantScalar):
                    # ConstantScalars can be added up immediately
                    new_value += operand.value
                else:
                    terms.append([1, [operand]])

            keyfunc = operator.itemgetter(1)

            terms.sort(key=keyfunc)

            # print("terms", terms)

            new_terms = []
            for key, group in itertools.groupby(terms, key=keyfunc):
                numeric_list, symbolic_list = list(zip(*group))

                numeric_sum = 0
                for number in numeric_list:
                    if isinstance(number, int):
                        # Test for int is ok because the only objects that are not
                        # ConstantScalars are the integers that were added a
                        # few lines above.
                        numeric_sum += number
                    else:
                        # number is ConstantScalar
                        numeric_sum += number.value

                # print("numeric_sum", numeric_sum)
                # print("symbolic_list", symbolic_list)

                # Constructing new terms.
                if numeric_sum == 1:
                    if len(symbolic_list[0]) > 1:
                        new_terms.append(ae.Times(symbolic_list[0]))
                    else:
                        new_terms.append(symbolic_list[0][0])
                elif numeric_sum == 0:
                    continue
                else:
                    new_terms.append(ae.Times(ae.ConstantScalar(numeric_sum), *symbolic_list[0]))

                # print("new_terms", new_terms)

            # The numeric value is only added if it is not zero, or if there are
            # no other terms.
            if new_value != 0:
                new_terms.append(ae.ConstantScalar(new_value))
            elif not new_terms:
                # new_values == 0, but new_terms is empty
                new_terms.append(ae.ConstantScalar(new_value))

            operands = new_terms 

            # print("after", expr)

        else:

            ############################
            # Dealing with non-scalar sums.

            terms = []
            for operand in expr.operands:
                if isinstance(operand, ae.Times):
                    scalars, non_scalars = operand.split_operands()
                    if scalars:
                        terms.append([scalars, non_scalars])
                    else:
                        terms.append([[1], non_scalars])
                else:
                    terms.append([[1], [operand]])

            keyfunc = operator.itemgetter(1)
            # def keyfunc(x):
            #     return x[1]

            terms.sort(key=keyfunc)

            # print("terms", terms)

            # Constructing new terms

            new_terms = []
            for key, group in itertools.groupby(terms, key=keyfunc):
                scalar_list, non_scalar_list = list(zip(*group))

                integers = []
                numeric_constants = []
                remaining_scalars = []

                for scalar in scalar_list:
                    if len(scalar) == 1:
                        if isinstance(scalar[0], int):
                            # It is ok to test for int because the only objects
                            # that are not ConstantScalars are the integers
                            # that were added a few lines above.
                            integers.append(scalar[0])
                        elif isinstance(scalar[0], ae.ConstantScalar):
                            # It is safe to assume that there are no products of
                            # multiple ConstantScalars because they were
                            # eliminated when the operands of the current were
                            # simplified. 
                            numeric_constants.append(scalar[0])
                        else:
                            # Single symbolic scalar.
                            remaining_scalars.append(scalar[0])
                    else:
                        # What remains are combinations of symbolic and
                        # numeric scalars. Nothing is done with them here.
                        # Instead, they are simplified when the new terms are
                        # created.
                        remaining_scalars.append(ae.Times(*scalar))

                # print("integers", integers)
                # print("numeric_constants", numeric_constants)
                # print("remaining", remaining_scalars)

                # print(integers)
                # Summing up all numeric values.
                numeric_sum = sum(integers)
                for numeric_constant in numeric_constants:
                    numeric_sum += numeric_constant.value

                # sum_of_scalars_list contains the sum of all scalars factors of
                # a certain expression. Example: for alpha A + 2 A, it is
                # [alpha, 2]

                if numeric_sum == 0:
                    if remaining_scalars:
                        sum_of_scalars_list = remaining_scalars
                    else:
                        # numeric scalars add up to zero and there are no
                        # symbolic ones, thus this term disappears
                        continue
                elif numeric_sum == 1 and not remaining_scalars:
                    # Numeric scalar is one and there are no remaining ones, so
                    # scalars are omitted entirely.
                    sum_of_scalars_list = []
                else:
                    sum_of_scalars_list = [ae.ConstantScalar(numeric_sum)] + remaining_scalars

                # print("plus", sum_of_scalars_list)        

                # Creating new terms.
                if sum_of_scalars_list:
                    if len(sum_of_scalars_list) > 1:
                        sum_of_scalars_expr = simplify(ae.Plus(*sum_of_scalars_list))
                        if isinstance(sum_of_scalars_expr, ae.ConstantScalar) and sum_of_scalars_expr.value == 0:
                            # same as: if sum_of_scalars_expr == ae.ConstantScalar(0):
                            # is it faster?
                            # After simplifying, it may happen that this sum
                            # ends up being zero. In that case, no new term is
                            # added.
                            continue
                        new_terms.append(ae.Times(sum_of_scalars_expr, *non_scalar_list[0]))
                    else:
                        new_terms.append(ae.Times(*(sum_of_scalars_list + non_scalar_list[0])))
                elif len(non_scalar_list[0]) > 1:
                    new_terms.append(ae.Times(*non_scalar_list[0]))
                else:
                    new_terms.append(non_scalar_list[0][0])
                
                # print("terms", new_terms)

            operands = new_terms

            ############################
            # Remove zero matrices in sums.

            remove = []
            for i, operand in enumerate(operands):
                if operand.has_property(properties.ZERO):
                    # print("operand", operand)
                    remove.append(i)

            for i in reversed(remove):
                del operands[i]


        return ae.Plus(*operands)
    else:
        # TODO anything missing?
        return expr

def invert(expr):
    """Returns the inverse of expr.

    IMPORTANT: This function modifies expr. Since it is possible that the
    outermost operator changes, it is possible that the variable originally
    passed to this function points to a wrong node in the expression tree.
    To avoid any problems, use this fuction as follows:

    expr = invert(expr)

    Remark:
    Times and Plus are commutative when applied to scalars. To ensure that
    testing equality of expressions always works, scalars are always sorted.
    An alternative might be to use _test_relation_commutative in
    Operator.__eq__.
    """
    if isinstance(expr, ae.Transpose):
        return invert_transpose(expr.operand)
    elif isinstance(expr, ae.Inverse):
        return expr.operand
    elif isinstance(expr, ae.InverseTranspose):
        return transpose(expr.operand)
    elif isinstance(expr, ae.Conjugate):
        return invert_conjugate(expr.operand)
    elif isinstance(expr, ae.ConjugateTranspose):
        return invert_conjugate_transpose(expr.operand)
    elif isinstance(expr, ae.InverseConjugate):
        return conjugate(expr.operand)
    elif isinstance(expr, ae.InverseConjugateTranspose):
        return conjugate_transpose(expr.operand)
    elif isinstance(expr, ae.Times):
        return _distribute_inverse(expr, ae.Inverse, invert, True)
    elif isinstance(expr, ae.Symbol):
        if isinstance(expr, ae.ConstantScalar):
            return ae.ConstantScalar(1/expr.value)
        elif expr.has_property(properties.ORTHOGONAL):
            # is it possible for an orthogonal matrix to be symmetric?
            return ae.Transpose(expr)
        elif expr.has_property(properties.UNITARY):
            return ae.ConjugateTranspose(expr)
        else:
            return ae.Inverse(expr)
    elif isinstance(expr, ae.Plus):
        return ae.Inverse(expr)
    elif isinstance(expr, matchpy.Wildcard): # For Wildcards
        return ae.Inverse(expr)
    elif isinstance(expr, ae.Equal):
        raise TransformationError("invert can not be applied to %r" % expr) 
    else:
        raise NotImplementedError("Error in expression.invert() while dealing with %r" % expr)

def transpose(expr):
    """Returns the transpose of expr.

    This function automatically distributes the transpose across operators.
    It does not change the canonical form of this expression.

    IMPORTANT: This function modifies expr. Since it is possible that the
    outermost operator changes, it is possible that the variable originally
    passed to this function points to a wrong node in the expression tree.
    To avoid any problems, use this fuction as follows:

    expr = transpose(expr)
    """
    if isinstance(expr, ae.Transpose):
        return expr.operand
    elif isinstance(expr, ae.Inverse):
        return invert_transpose(expr.operand)
    elif isinstance(expr, ae.InverseTranspose):
        return invert(expr.operand)
    elif isinstance(expr, ae.Conjugate):
        return conjugate_transpose(expr.operand)
    elif isinstance(expr, ae.ConjugateTranspose):
        return conjugate(expr.operand)
    elif isinstance(expr, ae.InverseConjugate):
        return invert_conjugate_transpose(expr.operand)
    elif isinstance(expr, ae.InverseConjugateTranspose):
        return invert_conjugate(expr.operand)
    elif isinstance(expr, ae.Times):
        scalars, non_scalars = expr.split_operands()
        operands = scalars + [transpose(operand) for operand in reversed(non_scalars)]
        return ae.Times(*operands)
    elif isinstance(expr, ae.Plus):
        operands = [transpose(operand) for operand in expr.operands]
        return ae.Plus(*operands)
    elif isinstance(expr, ae.Symbol):
        if expr.has_property(properties.SCALAR):
            return expr
        elif expr.has_property(properties.SYMMETRIC):
            return expr
        else:
            return ae.Transpose(expr)
    elif isinstance(expr, ae.ConstantScalar): # TODO remove this. ConstantScalar is a scalar symbol, so the symbol branch already does the right thing
        return expr
    elif isinstance(expr, matchpy.Wildcard): # For Wildcards
        return ae.Transpose(expr)
    else:
        raise TransformationError("transpose can not be applied to %r" % expr)

def invert_transpose(expr):
    """Returns the inverse-transpose of expr.

    IMPORTANT: Always use as expr = invert_transpose(expr).

    For some remarks, see docstring of invert().
    """
    if isinstance(expr, ae.Transpose):
        return invert(expr.operand)
    elif isinstance(expr, ae.Inverse):
        return transpose(expr.operand)
    elif isinstance(expr, ae.InverseTranspose):
        return expr.operand
    elif isinstance(expr, ae.Conjugate):
        return invert_conjugate_transpose(expr.operand)
    elif isinstance(expr, ae.ConjugateTranspose):
        return invert_conjugate(expr.operand)
    elif isinstance(expr, ae.InverseConjugate):
        return conjugate_transpose(expr.operand)
    elif isinstance(expr, ae.InverseConjugateTranspose):
        return conjugate(expr.operand)
    elif isinstance(expr, ae.Symbol):
        if expr.has_property(properties.SCALAR):
            return ae.Inverse(expr)
        elif expr.has_property(properties.SYMMETRIC):
            return ae.Inverse(expr)
        elif expr.has_property(properties.ORTHOGONAL):
            return expr
        elif expr.has_property(properties.UNITARY):
            return ae.Conjugate(expr)
        else:
            return ae.InverseTranspose(expr)
    elif isinstance(expr, ae.ConstantScalar):
        return ae.ConstantScalar(1/expr.value)
    elif isinstance(expr, ae.Plus):
        operands = [transpose(operand) for operand in expr.operands]
        return ae.Inverse(ae.Plus(*operands))
    elif isinstance(expr, ae.Times):
        return _distribute_inverse(expr, ae.InverseTranspose, invert_transpose, False)
    elif isinstance(expr, matchpy.Wildcard): # For Wildcards
        return ae.InverseTranspose(expr)
    else:
        raise TransformationError("invert_transpose can not be applied to %r" % expr)

def conjugate(expr):
    """Returns the conjugate of expr.

    IMPORTANT: Always use as expr = conjugate(expr).

    For some remarks, see docstring of invert().
    """
    if isinstance(expr, ae.Transpose):
        return conjugate_transpose(expr.operand)
    elif isinstance(expr, ae.Inverse):
        return invert_conjugate(expr.operand)
    elif isinstance(expr, ae.InverseTranspose):
        return invert_conjugate_transpose(expr.operand)
    elif isinstance(expr, ae.Conjugate):
        return expr.operand
    elif isinstance(expr, ae.ConjugateTranspose):
        return transpose(expr.operand)
    elif isinstance(expr, ae.InverseConjugate):
        return invert(expr.operand)
    elif isinstance(expr, ae.InverseConjugateTranspose):
        return invert_transpose(expr.operand)
    elif isinstance(expr, ae.Times):
        scalars, non_scalars = expr.split_operands()
        scalars = [conjugate(operand) for operand in scalars]
        scalars.sort()
        non_scalars = [conjugate(operand) for operand in non_scalars]
        operands = scalars + non_scalars
        return ae.Times(*operands)
    elif isinstance(expr, ae.Plus):
        operands = [conjugate(operand) for operand in expr.operands]
        return ae.Plus(*operands)
    elif isinstance(expr, ae.Symbol):
        # TODO rest for property real?
        # hermitian missing
        return ae.Conjugate(expr)
    elif isinstance(expr, ae.ConstantScalar):
        return ae.ConstantScalar(expr.value.conjugate())
    elif isinstance(expr, matchpy.Wildcard): # For Wildcards
        return ae.Conjugate(expr)
    else:
        raise TransformationError("conjugate can not be applied to %r" % expr)

def conjugate_transpose(expr):
    """Returns the conjugate-transpose of expr.

    IMPORTANT: Always use as expr = conjugate_transpose(expr).

    For some remarks, see docstring of invert().
    """
    if isinstance(expr, ae.Transpose):
        return conjugate(expr.operand)
    elif isinstance(expr, ae.Inverse):
        return invert_conjugate_transpose(expr.operand)
    elif isinstance(expr, ae.InverseTranspose):
        return invert_conjugate(expr.operand)
    elif isinstance(expr, ae.Conjugate):
        return transpose(expr.operand)
    elif isinstance(expr, ae.ConjugateTranspose):
        return expr.operand
    elif isinstance(expr, ae.InverseConjugate):
        return invert_transpose(expr.operand)
    elif isinstance(expr, ae.InverseConjugateTranspose):
        return invert(expr.operand)
    elif isinstance(expr, ae.Times):
        scalars, non_scalars = expr.split_operands()
        scalars = [conjugate(operand) for operand in scalars]
        scalars.sort()
        non_scalars = [conjugate_transpose(operand) for operand in reversed(non_scalars)]
        operands = scalars + non_scalars
        return ae.Times(*operands)
    elif isinstance(expr, ae.Plus):
        operands = [conjugate_transpose(operand) for operand in expr.operands]
        return ae.Plus(*operands)
    elif isinstance(expr, ae.Symbol):
        if expr.has_property(properties.SCALAR):
            return ae.Conjugate(expr)
        elif expr.has_property(properties.HERMITIAN):
            return expr
        else:
            return ae.ConjugateTranspose(expr)
    elif isinstance(expr, ae.ConstantScalar):
        return ae.ConstantScalar(expr.value.conjugate())
    elif isinstance(expr, matchpy.Wildcard): # For Wildcards
        return ae.ConjugateTranspose(expr)
    else:
        raise TransformationError("conjugate_transpose can not be applied to %r" % expr)    

def invert_conjugate(expr):
    """Returns the inverse-conjugate of expr.

    IMPORTANT: Always use as expr = invert_conjugate(expr).

    For some remarks, see docstring of invert().
    """
    if isinstance(expr, ae.Transpose):
        return invert_conjugate_transpose(expr.operand)
    elif isinstance(expr, ae.Inverse):
        return conjugate(expr.operand)
    elif isinstance(expr, ae.InverseTranspose):
        return conjugate_transpose(expr.operand)
    elif isinstance(expr, ae.Conjugate):
        return invert(expr.operand)
    elif isinstance(expr, ae.ConjugateTranspose):
        return invert_transpose(expr.operand)
    elif isinstance(expr, ae.InverseConjugate):
        return expr.operand
    elif isinstance(expr, ae.InverseConjugateTranspose):
        return transpose(expr.operand)
    elif isinstance(expr, ae.Symbol):
        if expr.has_property(properties.UNITARY):
            return ae.Transpose(expr)
        # TODO hermitian missing
        else:
            return ae.InverseConjugate(expr)
    elif isinstance(expr, ae.ConstantScalar):
        return ae.ConstantScalar(1/expr.value.conjugate())
    elif isinstance(expr, ae.Plus):
        operands = [conjugate(operand) for operand in expr.operands]
        return ae.Inverse(ae.Plus(*operands))
    elif isinstance(expr, ae.Times):
        return _distribute_inverse(expr, ae.InverseConjugate, invert_conjugate, True)
    elif isinstance(expr, matchpy.Wildcard): # For Wildcards
        return ae.InverseConjugate(expr)
    else:
        raise TransformationError("invert_conjugate can not be applied to %r" % expr) 
    
def invert_conjugate_transpose(expr):
    """Returns the inverse-conjugate-transpose of expr.

    IMPORTANT: Always use as expr = invert_conjugate_transpose(expr).

    For some remarks, see docstring of invert().
    """
    if isinstance(expr, ae.Transpose):
        return invert_conjugate(expr.operand)
    elif isinstance(expr, ae.Inverse):
        return conjugate_transpose(expr.operand)
    elif isinstance(expr, ae.InverseTranspose):
        return conjugate(expr.operand)
    elif isinstance(expr, ae.Conjugate):
        return invert_transpose(expr.operand)
    elif isinstance(expr, ae.ConjugateTranspose):
        return invert(expr.operand)
    elif isinstance(expr, ae.InverseConjugate):
        return transpose(expr.operand)
    elif isinstance(expr, ae.InverseConjugateTranspose):
        return expr.operand
    elif isinstance(expr, ae.Symbol):
        if expr.has_property(properties.SCALAR):
            return ae.InverseConjugate(expr)
        elif expr.has_property(properties.UNITARY):
            return expr
        elif expr.has_property(properties.HERMITIAN):
            return ae.Inverse(expr)
        else:
            return ae.InverseConjugateTranspose(expr)
    elif isinstance(expr, ae.ConstantScalar):
        return ae.ConstantScalar(1/expr.value.conjugate())
    elif isinstance(expr, ae.Plus):
        operands = [conjugate_transpose(operand) for operand in expr.operands] 
        return ae.Inverse(ae.Plus(*operands))
    elif isinstance(expr, ae.Times):
        return _distribute_inverse(expr, ae.InverseConjugateTranspose, invert_conjugate_transpose, False)
    elif isinstance(expr, matchpy.Wildcard): # For Wildcards
        return ae.InverseConjugateTranspose(expr)
    else:
        raise TransformationError("invert_conjugate_transpose can not be applied to %r" % expr) 

# @profile
def _distribute_inverse(expr, operator, operator_fct, reverse):
    """Distributes the inverse operators over a product.

    Additionally, scalars are moved to the front of the product.

    Distributing the inverse over products works similar for the operators 
    Inverse, InverseTranspose, InverseConjugate and InverseConjugateTranspose.
    This function helps to avoid code duplication.

    This function takes the following arguments:
    - expr: The expression over which the operator has to be distributed. It is
      assumed that is of type Times.
    - operator: The operator that has to be distributed. One of the following
      four: Inverse, InverseTranspose, InverseConjugate and
      InverseConjugateTranspose.
    - operator_fct: The corresponding function for the operator. One of the
      following four: invert, invert_transpose, invert_conjugate and
      invert_conjugate_transpose.
    - reverse: Boolean that indicates whether the order of the operand has to be
      revered or not. True for Inverse and InverseConjugate, False for
      InverseTranspose and InverseConjugateTranspose.
    """

    left = None
    right = None

    def split_operator(expr, operator):
        # Complex conjugation and transposition can always be distributed over
        # Times.
        if operator is ae.InverseTranspose:
            return ae.Inverse(transpose(expr))
        elif operator is ae.InverseConjugate:
            return ae.Inverse(conjugate(expr))
        elif operator is ae.InverseConjugateTranspose:
            return ae.Inverse(conjugate_transpose(expr))
        else:
            return operator(expr)

    scalars, non_scalars = expr.split_operands()
    scalars = [operator_fct(operand) for operand in scalars]
    scalars.sort()

    if len(scalars) == len(expr.operands):
        return ae.Times(*scalars)

    # print("distribute IN", expr)
    # print("SCALARS", scalars)
    # print(non_scalars)
    # print([op for op in non_scalars if op.has_property(properties.SQUARE)])

    # The operands are stored in a list. The operator can be distributed
    # over the expressions at the first level of the list. The nested list
    # contain those parts where the operator can not be distributed (because
    # the operands in those parts are not square).
    operands = []
    if expr.has_property(properties.SQUARE):
        # For square products, all square sub-sequences are collected in
        # sub-lists.
        # Example:
        # M1 = Matrix("M1", (8, 8))
        # M2 = Matrix("M2", (8, 16))
        # M3 = Matrix("M3", (16, 8))
        # M4 = Matrix("M4", (8, 8))
        # Times([M1, M2, M3, M4]) results in [M1, [M2, M3], M4]
        size = None # size is known. it has to be the size of the entire product. maybe this could be used to simplify things
        pos = None
        for n, operand in enumerate(non_scalars):
            if not size:
                if operand.has_property(properties.SQUARE):
                    operands.append(operand)
                else:
                    size = operand.size[0]
                    pos = n

            if size == operand.size[1]:
                operands.append(non_scalars[pos:n+1])
                size = None
    else:
        # Rectangular products are decomposed into
        # [<square factors>, [<rectangular part>], <square factors>]
        # The square parts can be empty.
        pos1 = None
        pos2 = None
        for n, operand in enumerate(non_scalars):
            if not operand.has_property(properties.SQUARE):
                pos1 = n
                break
        for n, operand in enumerate(reversed(non_scalars)):
            if not operand.has_property(properties.SQUARE):
                pos2 = len(non_scalars) - n
                break
        operands = non_scalars[0:pos1] + [non_scalars[pos1:pos2]] + non_scalars[pos2:]

    if reverse:
        operands.reverse()

    # The nested list representation is converted into an actual expression
    new_operands = []
    for operand in operands:
        if isinstance(operand, list):
            if len(operand) == 1:
                new_operands.append(operator_fct(operand[0]))
            else:
                new_operands.append(split_operator(ae.Times(*operand), operator))
        else:
            new_operands.append(operator_fct(operand))

    # new_operands = flatten_times(new_operands)
    # print("new_operands", new_operands)

    operands = scalars + new_operands
    # print("operands", operands)
    return ae.Times(*operands)

    # if not scalars and len(new_operands) == 1:
    #     expr = new_operands[0]
    # else:
    #     operands = scalars + new_operands
    #     expr = ae.Times(*operands)
    
    # print("distribute OUT", expr)
    # expr.to_dot_file()
    # return expr

def flatten_times(list):
    """Flattens Times in a list of expressions.

    Given a list [a, Times(b, c)], this function returns [a, b, c].
    """
    new_list = []
    for expr in list:
        if isinstance(expr, ae.Times):
            new_list.extend(expr.children)
        else:
            new_list.append(expr)
    return new_list