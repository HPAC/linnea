
from ..algebra.expression import Symbol, Transpose, Times, Equal, \
                                ConjugateTranspose, Inverse, \
                                InverseTranspose, Matrix, Plus

from ..algebra.transformations import transpose, invert, invert_transpose, simplify
from ..algebra.representations import to_SOP
from ..algebra.equations import Equations

from .gst import GST
from .utils import is_blocked

from collections import namedtuple, deque

from .graph.base.base import EdgeLabel

from .. import temporaries

import operator
import itertools
import copy
import math

import matchpy

def CSE_replacement_general(equations, expressions, expr_positions):

    if expressions:
        CSEs = find_CSEs_general(expressions, expr_positions)
    else:
        return []

    transformed_expressions = []

    # This is used to determine where to insert the equation for the
    # extracted CSE. It's placed directly before the first occurrence.
    min_eqn_idx = math.inf

    # This dictionary represents all possible combinations of different types.
    # Example: type_mapping[2][3] inverting (2) a transpose(inverse) (3) operator results
    # in a transposed operator (1).
    type_mapping = [[0, 1, 2, 3],
                    [1, 0, 3, 2],
                    [2, 3, 0, 1],
                    [3, 2, 1, 0]]

    # print(CSEs)

    tmp_type = None
    for CSE in CSEs:
        tmp = None
        eqn = None

        tmp_type, full_pos = CSE[0]
        eqn_idx, pos = full_pos

        # Creating temporary
        CSE_expr = equations[eqn_idx][pos]
        tmp = temporaries.create_tmp(CSE_expr, True)
        eqn = Equal(tmp, CSE_expr)

        # print(tmp)
        # print(eqn)

        equations_list = list(equations.equations)
        # print(equations_list)

        # Replacing all occurrences of this CSE
        replacements_per_equation = dict()
        for CSE_type, full_pos in CSE:
            eqn_idx, pos = full_pos

            min_eqn_idx = min(eqn_idx, min_eqn_idx)
            type = type_mapping[tmp_type][CSE_type]

            if type == 0:
                tmp_expr = tmp
            elif type == 1:
                tmp_expr = Transpose(tmp)
            elif type == 2:
                tmp_expr = Inverse(tmp)
            elif type == 3:
                tmp_expr = InverseTranspose(tmp)

            replacements_per_equation.setdefault(eqn_idx, []).append((pos, tmp_expr))

        for eqn_idx, replacements in replacements_per_equation.items():
            equations_list[eqn_idx] = matchpy.replace_many(equations_list[eqn_idx], replacements)

        # Inserting new equation for extracted CSE.
        # It is inserted right before the first occurrence of the CSE.
        equations_list.insert(min_eqn_idx, eqn)
        new_equations = Equations(*equations_list)
        new_equations = new_equations.to_normalform()
        # new_equations._cleanup()
        # print(new_equations)
        transformed_expressions.append((new_equations, new_equations.metric(), EdgeLabel()))

    return transformed_expressions


def CSE_replacement_times(equations, products, product_positions):

    # print("#####")
    # print(products)
    # print(eqn_indices)
    if products:
        CSEs = find_CSEs_times(products)
    else:
        return []
    # print("CSEs", CSEs)
    # print(product_positions)

    transformed_expressions = []

    # This dictionary represents all possible combinations of different types.
    # Example: type_mapping[2][3] inverting (2) a transpose(inverse) (3) operator results
    # in a transposed operator (1).
    type_mapping = [[0, 1, 2, 3],
                    [1, 0, 3, 2],
                    [2, 3, 0, 1],
                    [3, 2, 1, 0]]
    
    # Stores the type of the operands of tmp.
    tmp_type = None

    for positions, length in CSEs:
        # print("##")
        # print(positions)
        # Create tmp
        seq_idx, occurrences = next(iter(positions.items()))
        CSE_pos = occurrences[0]
        # print(CSE_pos)
        operands = products[seq_idx].operands[CSE_pos.pos:CSE_pos.pos+length]    
        CSE_expr = Times(*operands)
        if is_blocked(CSE_expr):
            continue
        tmp_type = CSE_pos.type
        if tmp_type == 1:
            CSE_expr = transpose(CSE_expr)
        elif tmp_type == 2:
            CSE_expr = invert(CSE_expr)
        elif tmp_type == 3:
            CSE_expr = invert_transpose(CSE_expr)

        tmp = temporaries.create_tmp(CSE_expr, True)
        eqn = Equal(tmp, CSE_expr)

        replacements_per_equation = dict()
        equations_list = list(equations.equations)
        min_seq_idx = math.inf
        for seq_idx, occurrences in positions.items():
            CSE_eqn_idx, CSE_path = product_positions[seq_idx]
            # print(product_positions[seq_idx])
            operands = equations_list[CSE_eqn_idx][CSE_path].operands.copy()

            # print("occurrences", occurrences)
            # Remember: occurrences are sorted by pos.
            # print(equations_list)
            for CSE_pos in reversed(occurrences):
                min_seq_idx = min(seq_idx, min_seq_idx)
                
                # "Empty" replacements are added to remove operands.
                for i in range(CSE_pos.pos+1, CSE_pos.pos+length):
                    # print(CSE_path+(i,))
                    replacements_per_equation.setdefault(CSE_eqn_idx, []).append((CSE_path+(i,), []))

                del operands[CSE_pos.pos:CSE_pos.pos+length]
                type = CSE_pos.type
                if type == 0:
                    tmp_expr = tmp
                elif type == 1:
                    tmp_expr = Transpose(tmp)
                elif type == 2:
                    tmp_expr = Inverse(tmp)
                elif type == 3:
                    tmp_expr = InverseTranspose(tmp)

                replacements_per_equation.setdefault(CSE_eqn_idx, []).append((CSE_path + (CSE_pos.pos,), tmp_expr))
                operands.insert(CSE_pos.pos, tmp_expr)
        
        # print(replacements_per_equation)
        for eqn_idx, _replacements in replacements_per_equation.items():
            equations_list[eqn_idx] = matchpy.replace_many(equations_list[eqn_idx], _replacements)

        # Inserting new equation for extracted CSE.
        # It is inserted right before the first occurrence of the CSE.
        equations_list.insert(product_positions[min_seq_idx][0], eqn)
        new_equations = Equations(*equations_list)
        new_equations = new_equations.to_normalform()

        transformed_expressions.append((new_equations, new_equations.metric(), EdgeLabel()))
    return transformed_expressions

def CSE_replacement_plus(equations, sums, sum_positions):

    if sums:
        CSEs = find_CSEs_plus(sums)
    else:
        return []

    # print("CSEs", CSEs)

    transformed_expressions = []

    # print(equations)
    for tmp_description, remove, tmp_type in CSEs:
        # print(tmp_description, remove, tmp_type)
        min_seq_idx = math.inf

        operand_lists = [sum.operands.copy() for sum in sums]
        equations_list = list(equations.equations)

        seq_idx, positions = tmp_description
        tmp_operands = [operand_lists[seq_idx][pos] for pos in positions]
        CSE_expr = Plus(*tmp_operands)
        # print(CSE_expr)
        tmp = temporaries.create_tmp(CSE_expr, True)
        eqn = Equal(tmp, CSE_expr)
        # print(tmp.name)
        # print(tmp.get_equivalent())
        # print(eqn)

        for seq_idx, positions in remove.items():
            min_seq_idx = min(seq_idx, min_seq_idx)
            positions.sort(reverse=True)
            # print(seq_idx, positions)
            for pos in positions:
                # print(sums_copy[seq_idx])
                # print(sums_copy[seq_idx].operands, pos)
                operand_lists[seq_idx].pop(pos)

        for seq_idx, types in tmp_type.items():
            # tmp_expr = None
            for type in types:
                tmp_expr = None
                if type == 0:
                    tmp_expr = tmp
                elif type == 1:
                    tmp_expr = Transpose(tmp)

                operand_lists[seq_idx].append(tmp_expr)

        for seq_idx, operands in enumerate(operand_lists):
            CSE_eqn_idx, CSE_path = sum_positions[seq_idx]
            # TODO use replace_many
            equations_list[CSE_eqn_idx] = matchpy.replace(equations_list[CSE_eqn_idx], CSE_path, Plus(*operands))

        # for sum in sums_copy:
        #     print(sum.operands)
        # print(repr(equations_list[1]))
        equations_list.insert(sum_positions[min_seq_idx][0], eqn)

        new_equations = Equations(*equations_list)
        new_equations = new_equations.to_normalform()

        transformed_expressions.append((new_equations, new_equations.metric(), EdgeLabel()))
    return transformed_expressions


def _is_prefix(p1, p2):
    """Tests if one position is a prefix of another.

    Those positions are expected two be tuples of two values:
    - An integer, the equation index.
    - A list of integers, the position of a subexpression in an expression.
    """
    if p1[0] != p2[0]:
        return False
    else:
        return all(x == y for x, y in zip(p1[1], p2[1]))

def find_CSEs_general(expressions, expr_positions):
    """Finds general common subexpressions.

    As input, this function requires:
    - A list of expressions.
    - A list containing the positions of those expressions. It is assumed that
      the expressions and their corresponding positions are in the same postions
      in those lists.

    This function returns a list of lists. Each list represents one CSE. Those 
    lists contain tuples. Each tuple represents one occurrence of this CSE:
    - The first entry is an integer, specifying the type of this occurrence of
      the CSE.

      0 unmodified
      1 transposed
      2 inverted
      3 inverted and transposed
    - The second entry is the position of this occurrence of the CSE. It is
      taken, unmodified, from the input argument expr_positions.

    The idea for this function is the following: In principle, the generalized
    suffix tree could be used. However, here, we are not dealing with sequences
    of expressions, but with single expressions, so the GST adds a lot of
    unnecessary overhead. Thus, the following is done, inspired by what the
    GST does if all sequences have length one.
    An initially empty list of expressions is created (expr_list), together with
    a second list for the positions (pos_list). For each new expression, it is
    first checked if an equal expression is already in the list. If this is not
    the case, the new expression is appended to the list of expressions, and its
    position is appended to the list of positions (inside a list of length one).

    If the expression is already in the list, the new position is appended to
    the list of positions for that expression. This results in something like
    this:
    expr_list = [expr1, expr2]
    pos_list = [[p1], [p2, p3]]
    Here, expr2 is a common supexpression because it appears two times (in
    positions p2 and p3).
    """

    # expressions_T = [transpose(copy.deepcopy(expr)) for expr in expressions]
    # expressions_INV = [invert(copy.deepcopy(expr)) for expr in expressions]

    expressions_T = [transpose(expr) for expr in expressions]
    expressions_INV = [invert(expr) for expr in expressions]

    expr_list = []
    pos_list = []


    # print(expressions)
    for type, exprs in enumerate([expressions, expressions_T, expressions_INV]):
        for n1, expr1 in enumerate(exprs):
            for n2, expr2 in enumerate(expr_list):
                if expr1 == expr2:
                    expr_pos = expr_positions[n1]
                    # With is_prefix, two things are checked:
                    # First, it is tested whether this positions is already
                    # in the list of positions. This can happen if an
                    # expression is symmetric (because then, expr = expr^T).
                    # Second, it is tested whether one positions is a prefix of
                    # another expression. This can happen for an expression like
                    # Inverse(Plus(A, B)), where
                    # invert(Inverse(Plus(A, B))) == Plus(A, B)
                    if not any(_is_prefix(pos, expr_pos) for _, pos in pos_list[n2]):
                        pos_list[n2].append((type, expr_pos))
                        break
            else:
                expr_list.append(expr1)
                pos_list.append([(type, expr_positions[n1])])

    remove = []
    for n, potential_CSE in enumerate(pos_list):
        if len(potential_CSE) > 1:
            type0 = potential_CSE[0][0]
            # If all positions of a potential CSE have the same type!=0, then
            # this is also a CSE where all postions have type == 0. To avoid
            # duplicates, all potential CSEs where all positions have the same
            # type!=0 are removed.
            if type0 == 0:
                continue
            for type, positions in potential_CSE[1:]:
                if type != type0:
                    break
            else:
                remove.append(n)
        else:
            # Expressions that only show up once can not be common subexpressions.
            remove.append(n)


    for n in reversed(remove):
        # del expr_list[n]
        del pos_list[n]

    # for p, e in zip(pos_list, expr_list):
    #     print(p, e)
    # print("\n")

    return pos_list

CSEPosition = namedtuple("CSEPosition", ["pos", "type"])

# @profile
def find_CSEs_times(products):
    """Finds all replaceable CSEs in the input (products).

    The input has to be a list of products, i.e. Times objects.

    This function returns a list of tuples. Each tuple contains
    - a dictionary
    - the length of this CSE.

    The dictionary maps integers to lists of named tuples CSEPosition.
    Example:
    {0: [CSEPosition(pos=0, type=0), CSEPosition(pos=2, type=1)]}
    The keys of the dictionary represent positions in the input list (products).
    Thus, this means that the following occurrences of this CSE were found in
    the first item of the input (the first product).

    Each CSEPosition object describes one occurrence. pos is the position in the
    list of operands of that product (the list of operands, starting at 0). type
    specifies if this occurrence is
    0 unmodified
    1 transposed
    2 inverted

    The list in the dictionary are sorted by pos.
    """
    n = len(products)
    counter = 0

    # Here, all operands that are not simple are replaced with a placeholder.
    # This is necessary because otherwise, non-simple expression can cause
    # problems. Example: Times(A, B, Inv(Times(A, B))). If this is inverted,
    # Times(A, B) basically moves into the outer times and messes up all
    # indices.
    _products = []
    for _product in products:
        new_operands = None
        for _n, _operand in enumerate(_product.operands):
            if not _is_simple_times_CSE(_operand):
                if not new_operands:
                    new_operands = _product.operands.copy()
                # TODO if this ever causes problems, give this matrix a random
                # name (or maybe just str(_operand.size)). It could cause
                # problems because multiple matrices with the same name are
                # created.
                # I'm not sure if using the actual size is necessary.
                name = "".join(["dummy", str(counter)])
                counter += 1
                new_operands[_n] = Matrix(name, _operand.size)
                # _product.operands[_n] = Matrix(name, _operand.size)
        if new_operands:
            _products.append(Times(*new_operands))
        else:    
            _products.append(_product)

    products = _products

    # print(products)
    products_T = [transpose(product) for product in products]
    # print(products_T)

    # Mathematically, this is not the correct way to compute the inverse of a
    # general product. However, for the detection of common subexpression to
    # work correctly, the order of the factors has to be reversed, and the
    # inverse has to be applied to each factor individually.
    # (The mathematically correct way to do this would be to only move those
    # factors out of the inverse that are invertible.)
    products_INV = [Times(*(invert(factor) for factor in reversed(product.operands))) for product in products]

    gst = GST()
    for product in products + products_T + products_INV:
        if isinstance(product, Symbol):
            gst.add_sequence([product])
        else:    
            gst.add_sequence(product.operands)

    # print(products)
    # print(products_T)
    # print(products_INV)
    CSEs = gst.find_all_CSEs()
    # gst.to_dot_file()
    # quit()
    # print(CSEs)
    actual_CSEs = []

    # First, we transform the CSEs in different versions (i.e. transposed,
    # inverted) of the sequence into CSEs in the original, unmodified sequence,
    # with additional information about their type (i.e. unmodified, transposed,
    # …).
    for positions, length in CSEs:
        # Computing type based on seq_idx.
        # type = floor(seq_idx/n)
        # type = 0 : unmodified occurrence of CSE
        # type = 1 : transposed occurrence of CSE
        # type = 2 : inverted occurrence of CSE
        types = [seq_idx//n for seq_idx in positions.keys()]
        if all(type == 0 for type in types):
            # If all occurences appear in unmodified sequences, there is not much
            # to do.
            actual_positions = dict()
            for seq_idx in positions:
                for pos in positions[seq_idx]:
                    actual_positions.setdefault(seq_idx, []).append(CSEPosition(pos, 0))
            actual_CSEs.append((actual_positions, length))

        # Remark: If all occurences appear in the same modified sequence, we
        # ignore that CSE because it's also a CSE in the unmodified sequence.
        # This is the case if all entries in types are the same, but unequal 0.
        # Example: A B C A B           CSE: A B
        #          B^T A^T C^T B^T A^T CSE: B^T A^T

        for type in types[1:]:
            if types[0] != type:
                # Here, we look at CSEs that appear in different version of the
                # sequence. In this case, we actually have to map the
                # information to the unmodified sequence, recalculating the 
                # actual position.
                actual_positions = dict()
                for seq_idx in positions:
                    for pos in positions[seq_idx]:
                        _type = seq_idx // n
                        actual_seq_idx = seq_idx % n
                        if _type > 0:
                            actual_pos = len(products[actual_seq_idx].operands) - length - pos
                            if actual_pos < 0:
                                print(positions, length)
                        else:
                            actual_pos = pos
                        actual_positions.setdefault(actual_seq_idx, []).append(CSEPosition(actual_pos, _type))
                actual_CSEs.append((actual_positions, length))
                break

    # for CSE in actual_CSEs:
    #     print(CSE, is_replaceable_CSE(CSE))
    # print("actual_CSEs", actual_CSEs)

    # Returns only those CSEs that are replaceable.
    # Careful: is_replaceable_CSE intentionally sorts the lists of occurrences.
    actual_CSEs = list(filter(is_replaceable_CSE, actual_CSEs))

    # print("actual_CSEs 1", actual_CSEs)

    # Removing duplicates
    # Here, duplicates means that all positions are the same. Types might be
    # different (as far as I know, the algorithm never generates true
    # duplicates).
    # remove = []
    # for cse1, cse2 in itertools.combinations(actual_CSEs, 2):
    #     if equal_CSE_times(cse1, cse2):
    #         remove.append(cse1)

    # for cse in remove:
    #     actual_CSEs.remove(cse)

    # print("actual_CSEs", actual_CSEs)

    # Removing CSEs that are not maximal, i.e. they can be extended, and the
    # extended expression is still a CSE with the same number of occurences.
    remove = []
    for cse1, cse2 in itertools.product(actual_CSEs, repeat=2):
        if cse1 == cse2:
            continue
        elif not_max_CSE_times(cse1, cse2) or sub_CSE_times(cse1, cse2):
            remove.append(cse1)

    for cse in remove:
        # In some cases, certain CSEs are detected as "not maximal" multiple
        # times (this relation is transitive).
        try:
            actual_CSEs.remove(cse)
        except ValueError:
            pass

    # print("actual_CSEs", actual_CSEs)

    return actual_CSEs

def equal_CSE_times(cse1, cse2):
    """Test if two CSEs are equal.

    This function can NOT be replaced by the == operator. This function does
    not consider the type for equality. CSEs that just differ in the type
    can lead to exactly the same replacement at the moment (this has
    something to do with what type is chosen as the "reference type" for
    the extracted CSE. Just one type is used, not all possible types.)
    """
    positions1, length1 = cse1
    positions2, length2 = cse2
    if length1 != length2:
        return False
    elif len(positions1) != len(positions2):
        return False
    else:
        for seq_idx, occurrences1 in positions1.items():
            try:
                occurrences2 = positions2[seq_idx]
            except KeyError:
                return False
            else:
                if len(occurrences1) != len(occurrences2):
                    return False
                elif any(CSE_pos1.pos != CSE_pos2.pos for CSE_pos1, CSE_pos2 in zip(occurrences1, occurrences2)):
                    return False
        return True

def sub_CSE_times(cse1, cse2):
    """Tests if cse1 is a true "sub-CSE" of cse2.

    cse1 is a "sub-CSE" if all occurrences of cse1 are also occurrences in
    cse2 (not considering type). This function returns False if both CSEs
    are equal.

    Intuitively, both CSEs are the same expression, but cse1 is missing some
    occurrences.
    """
    positions1, length1 = cse1
    positions2, length2 = cse2
    if length1 != length2:
        return False
    else:
        # There has to be at least one s1 that is a subset of s2,
        # otherwise, both CSEs might be equal.
        subset_exists = False
        for seq_idx, occurrences2 in positions2.items():
            try:
                occurrences1 = positions1[seq_idx]
            except KeyError:
                # If there are no occurrences for that seq_idx, this
                # means that the set of occurrences is empty. The empty set
                # is subset of every other set.
                subset_exists = True
            else:
                if len(occurrences1) > len(occurrences2):
                    return False
                else:
                    s1 = set(CSE_pos1.pos for CSE_pos1 in occurrences1)
                    s2 = set(CSE_pos2.pos for CSE_pos2 in occurrences2)
                    # print("s1", s1)
                    # print("s2", s2)
                    # If s1 is not a subset of s2, cse1 can certainly not be a
                    # subset of cse2.
                    # Keep in mind that this is not the same as s1 > s2
                    if not s1 <= s2:
                        return False
                    if s1 < s2:
                        subset_exists = True
        return subset_exists

# def sub_CSE_times(cse1, cse2):
#     """Tests if cse1 is a "sub-CSE" of cse2.

#     cse1 is a "sub-CSE" if all occurrences of cse1 are also occurrences in
#     cse2 (not considering type).

#     Intuitively, both CSEs are the same expression, but cse1 is missing some
#     occurrences.
#     """
#     positions1, length1 = cse1
#     positions2, length2 = cse2
#     if length1 != length2:
#         return False
#     else:
#         for seq_idx, occurrences1 in positions1.items():
#             try:
#                 occurrences2 = positions2[seq_idx]
#             except KeyError:
#                 return False
#             else:
#                 if len(occurrences1) > len(occurrences2):
#                     return False
#                 else:
#                     s1 = set(CSE_pos1.pos for CSE_pos1 in occurrences1)
#                     s2 = set(CSE_pos2.pos for CSE_pos2 in occurrences2)
#                     if s1 > s2:
#                         return False
#         return True



def not_max_CSE_times(cse1, cse2):
    """Tests if cse1 is a not maximal with respect to cse2.

    Maximal means that if a CSE is extended, then one or more occurrences
    are lost. See "Structural Analysis of Gapped Motifs of a String" by Esko
    Ukkonen, section "Representation by Maximal Non–gapped Motifs"

    Here, we simply test this by pairwise comparison. If cse2 is "larger"
    than cse1, than cse1 can not be maximal.

    The idea is to test if two CSEs which have different length have the
    same positions modulo a small offset (the offset can not be larger than
    the difference between the length).
    """
    
    positions1, length1 = cse1
    positions2, length2 = cse2
    
    if length1 < length2 and len(positions1) == len(positions2):
        length_difference = length2 - length1
        for offset in range(length_difference + 1):
            match = True
            for key1 in positions1.keys():
                try:
                    p2 = positions2[key1]
                except KeyError:
                    # In this case, cse1 is certainly not a subset of cse2
                    return False

                p1 = positions1[key1]
                if len(p1) != len(p2):
                    return False

                for occurrence1, occurrence2 in zip(p1, p2):
                    pos2 = None
                    if occurrence2.type == 0:
                        pos2 = occurrence2.pos + offset
                    else:
                        pos2 = occurrence2.pos + length_difference - offset
                    if not (occurrence1.type == occurrence2.type and occurrence1.pos == pos2):
                        match = False
            if match:
                return True
    return False

def _sequence_intersection(sequences):
    """Computes all intersections of all sequences.

    sequences has to be a list of sorted lists. This functions computes all
    intersections of those lists. 

    The idea is to simultaneously traverse all lists, collecting elements that
    appear in multiple lists. This is done by maintaining a list of pointers,
    where each pointer points to one element of a different list (the current
    elements). Initially, all pointers point to the first element.

    1. The smallest element of the current elements is determined.
    2. Those elements that are equal to this smallest element are "collected".
    3. The pointers pointing to those elements are moved forward by one.
    4. This process is repeated, starting at 1., until all lists are completely
       traversed.

    The intersections are returned as a dictionary.
    - The keys of this dictionary are tuples of integers. Those integers are the
      indices in sequences that are part of this intersection. Example: (1, 2)
      means that sequences[1] and sequences[2] have some elements in common.
    - The values of this dictionary are lists of lists. Each inner list is one
      set of indices, where those positions that also show up in the key are the
      indices of that element that is in the intersection of those list.

    Example:
    Input: [["a", "b", "c"],["b", "c", "d"],[1, "c", "d"]]
    Output: {(0, 1): [[1, 0, 1]], (0, 1, 2): [[2, 1, 1]], (1, 2): [[2, 2, 2]], (0, 2): [[0, 0, 0]]}

    (0, 1): [[1, 0, 1]] means the following: From (0, 1) it follows that the
    first two sequences have an element in common (it is "b"). Thus, we only
    care about position 0 and 1 in l = [1, 0, 1]. l[0]=1 tells us that this
    element of this intersection has the index 1 in sequences[0]. It is "b". The
    other occurrence is l[1]=0, that is, the first position in sequences[1].
    """

    n = len(sequences)

    pointer = [0 for _ in range(n)]
    not_end = [True for _ in range(n)]

    CSEs = dict()

    while True:

        # print("P:", pointer)

        # Finding the index k of the first sequences that has some unprocessed
        # elements left.
        k = not_end.index(True)
        min_val = sequences[k][pointer[k]]
        min_sequences = [k]

        # It is sufficient to start at k+1 because for all i < k, not_end[i] is
        # False.
        for i in range(k+1, n):
            if not_end[i]:
                val = sequences[i][pointer[i]]
                # print(val, min_val)
                if val < min_val:
                    # print("smaller")
                    min_sequences = [i]
                    min_val = val
                elif val == min_val:
                    # print("equal")
                    min_sequences.append(i)

        # print(min_sequences)

        if len(min_sequences) > 1:
            CSEs.setdefault(tuple(min_sequences), []).append(tuple(copy.deepcopy(pointer)))
            # CSEs.setdefault(tuple(min_sequences), []).append((min_val, copy.deepcopy(pointer)))

            # This loop generates additional subsets. Example: If sequences 1, 2
            # and 3 have an intersection, then (1, 2), (1, 3) and (2, 3) have
            # intersections as well. They are generated here.
            for i in range(2, len(min_sequences)):
                for min_sequences_subset in itertools.combinations(min_sequences, i):
                    CSEs.setdefault(tuple(min_sequences_subset), []).append(tuple(copy.deepcopy(pointer)))

        # Advancing pointers pointing to minimal elements.
        for i in min_sequences:
            if pointer[i] + 1 == len(sequences[i]):
                not_end[i] = False
            else:
                pointer[i] += 1

        # print(not_end)
        if not any(not_end):
            break

        
        # print(min_val)
        # print(min_sequences)

    return CSEs

class GraphNode(object):
    """Node for the graph that represents the relation between intersections.

    """
    def __init__(self, seq_indices, positions):
        self.seq_indices = seq_indices # frozenset
        self.positions = positions # list of lists
        self.supersets = []
        self.subsets = []

def _find_CSE_nodes(node):
    """

    This function traverses the subset graph starting at "node", move to all
    subsets. Whenever a subset is found with more than one position, it stops
    and returns that node as a CSE node.

    This way, by traversing the graph from top to bottom and stopping as soon
    as a CSE node is found, we don't detect CSEs where all occurrences are a
    subsets of the occurrences of another CSE (this is the relation that this
    graph represents).
    """
    if len(node.positions) > 1:
        yield node
    else:
        for n in node.subsets:
            yield from _find_CSE_nodes(n)

def _sort(l):
    """Returns a sorted copy of l plus a permutation list.

    The permutation list describes how the original list was reordered during
    sorting.
    """ 
    return zip(*sorted(enumerate(l), key=operator.itemgetter(1)))

def find_CSEs_plus(sums):
    n = len(sums)

    tmp_sequences = [copy.copy(sum.operands) for sum in sums]
    # Careful! transpose() sorts. Here however, I have to keep track of how the
    # summands change positions, so I have to apply transpose in a way that it
    # doesn't sort.
    tmp_sequences.extend([[transpose(copy.deepcopy(summand)) for summand in sum.operands] for sum in sums])

    pos_map = dict()
    sequences = []
    for seq_idx, sequence in enumerate(tmp_sequences):
        # Sequences are sorted to allow efficient computation of the
        # intersections. The pos_map stores the original positions of all
        # expressions.
        pos_map[seq_idx], sequence = _sort(sequence)
        sequences.append(sequence)

    # print("sequences", sequences)

    intersections = _sequence_intersection(sequences)
    # print(pos_map)
    # print("intersections", intersections)

    nodes = []
    for seq_indices, positions in intersections.items():

        # Here, we are dealing with the ugly special case. Notice that this
        # still doesn't work in all cases.
        # This is one example of the ugly special case
        # Example: A + B + A^T + B^T + C
        # The problem is that the intersection of that sum and its
        # transpose contains the CSE A + B (or A^T + B^T, or A + B^T ...)
        # twice (this does not happen if there is one additional occurrence
        # of this CSE in some other sum. However, it potentially gets even
        # worse when there are two sums like that. Note that A + B can not
        # show up twice in one sum because it will be simplified to
        # 2 A + 2 B).
        # Solution (at least for the case where there is one such sum):
        # We keep track of how many occurences there are in each orignal
        # sequence. sequence_hit_count stores how often each sequence is "hit"
        # If all sequences are either hit twice (or none), then this is a
        # special case. For some of those cases, the problem can be solved by 
        # discarding the second half of the positions. This works because the
        # expressions were sorted, so the order is always [transposed, not transposed]
        sequence_hit_count = [0 for _ in range(n)]

        for seq_idx in seq_indices:
            sequence_hit_count[seq_idx % n] += 1

        # print("sequence_hit_count", seq_indices, sequence_hit_count, positions)

        if all((hits == 2 or hits == 0) for hits in sequence_hit_count):
            positions = positions[:len(positions)//2]
            pass

        # print("sequence_hit_count", seq_indices, sequence_hit_count, positions)
        nodes.append(GraphNode(frozenset(seq_indices), positions))

    # TODO this graph is similar a lattice https://en.wikipedia.org/wiki/Lattice_(order)

    # Here, we are constructing the graph that respresents the subset/superset
    # relation between the intersections. That is, the relation between the
    # sequence indices.
    for node1, node2 in itertools.combinations(nodes, 2):
        if node1.seq_indices < node2.seq_indices:
            node1.supersets.append(node2)
            node2.subsets.append(node1)
        elif node1.seq_indices > node2.seq_indices:
            node1.subsets.append(node2)
            node2.supersets.append(node1)

    # A dictionary is used to make sure that no node is added more than ones.
    # Nodes can be uniquely identified by seq_indices.
    CSE_nodes = dict()

    for node in nodes:
        if not node.supersets:
            for nd in _find_CSE_nodes(node):
                seq_indices = list(nd.seq_indices)
                type = seq_indices[0] // n
                # If all occurences have the same type != 0, we also detect a CSE where
                # all occurrences have type 0. Thus, those CSEs can be ignored.
                # TODO It should be possible to discard those cases before contructing
                # the graph.
                if type != 0 and all(idx //n == type for idx in seq_indices[1:]):
                    continue
                else:
                    # print("HERE")
                    # Node is only added if its not in CSE_nodes.
                    CSE_nodes.setdefault(nd.seq_indices, nd)

    type_mapping = [[0, 1],
                    [1, 0]]

    # print(CSE_nodes)
    # print(pos_map)

    # Constructing the output representation fo the CSEs.
    CSEs = []
    for node in CSE_nodes.values():
        positions = node.positions

        # print("NODE:", node.seq_indices)   

        # We don't care about CSEs with less than two operands.
        if len(positions) > 1:
            remove = dict()
            tmp_types = dict()
            tmp = None
            reference_type = None
            for seq_idx in node.seq_indices:
                type = seq_idx // n
                actual_seq_idx = seq_idx % n
                actual_positions = [pos_map[seq_idx][position[seq_idx]] for position in positions]
                if not tmp:
                    reference_type = type
                    tmp = (actual_seq_idx, actual_positions)

                tmp_types.setdefault(actual_seq_idx, []).append(type_mapping[type][reference_type])
                
                remove.setdefault(actual_seq_idx, []).extend(actual_positions)
            # print(node.seq_indices)
            # print(positions)
            # print("tmp:", tmp)
            # print("Remove:", remove)
            # print("tmp_types", tmp_types)
            # if there are duplicates in remove, something is wrong
            duplicates = False
            for _, _positions in remove.items():
                seen = set()
                for _pos in _positions:
                    if _pos in seen:
                        duplicates = True
                        break
                    else:
                        seen.add(_pos)
                if duplicates:
                    break
            if not duplicates:
                CSEs.append((tmp, remove, tmp_types))

    return CSEs


def is_replaceable_CSE(CSE):
    """Test if a CSE is replaceable.

    Replaceable means that there are no occurrences that overlap with itself.
    Example: Sequence AAAA
             CSE 1    AAA
             CSE 2     AAA
    Clearly, occurrences of AAA  overlap with itself, so replacing both is not
    possible.
    Of course, if there is a third occurrence somewhere else, it would be
    possible to replace CSEs 1 and 3 or 2 and 3. We do not take this into
    account because we expect that cases like that are rare.
    In general, a subset of all CSEs might be replaceable. In order to explore
    every possibility, one had to inspect all subsets. It would also be possible
    to construct just those subsets that are replaceable, but the question is if
    that would be worth the effort.
    """
    positions, length = CSE
    for seq_idx, occurrences in positions.items():
        if len(occurrences) > 1:
            occurrences.sort(key=lambda x: x.pos)
            for n, CSE_pos in enumerate(occurrences[:-1]):
                if CSE_pos.pos + length > occurrences[n+1].pos:
                    return False
    return True

def is_simple_times_CSE(expr):
    """Test if expr is sufficiently simple for CSE replacement."""
    return len(expr.operands)>1 and all(_is_simple_times_CSE(operand) for operand in expr.operands)

def _is_simple_times_CSE(expr):
    # TODO what is missing?
    if isinstance(expr, Symbol) \
    or (isinstance(expr, Transpose) and isinstance(expr.operand, Symbol)) \
    or (isinstance(expr, ConjugateTranspose) and isinstance(expr.operand, Symbol)) \
    or (isinstance(expr, Inverse) and isinstance(expr.operand, Symbol)) \
    or (isinstance(expr, InverseTranspose) and isinstance(expr.operand, Symbol)):
        return True
    else:
        return False
