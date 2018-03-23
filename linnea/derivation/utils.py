
import matchpy
import math

def apply_kernel_with_context(expr, many_to_one_matcher):
    """Applies one kernel to an expressions.

    This function applies the first matching kernel of the many_to_one_matcher
    to expr and returns the replacement, as well as the matched kernel.

    The functions only searches for matches at the root of expr, and it expects
    many_to_one_matcher to use the patterns with context variables
    (Kernel.pattern_with_context).

    Args:
        expr (Expression): The expression where the kernel should be applied.
        many_to_one_matcher (ManyToOneMatcher): Matcher for kernel patterns with
            context. The label has to be the kernel.

    Returns:
        A tuple containing
        - the replacement (Expression)
        - the matched kernel (MatchedKernel)
        If no match was found, the replacement is the input expression and
        instead of the matched kernel, None is returned.
    """

    for kernel, substitution in many_to_one_matcher.match(expr):
        matched_kernel = kernel.set_match(substitution, True)
        expr = matched_kernel.replacement
        return expr, matched_kernel
    return (expr, None)


def apply_kernel_anywhere(expr, discrimination_net):
    """Applies one kernel to an expressions.

    This function applies the first matching kernel of the discrimination_net to
    expr and returns the replacement, as well as the matched kernel.

    The functions searches for matches anywhere in expr, and it expects
    discrimination_net to use the pattern without context variables
    (Kernel.pattern).

    Args:
        expr (Expression): The expression where the kernel should be applied.
        discrimination_net (DiscriminationNet): Discrimination net for kernel
        patterns without context. The final_label has to be the kernel.

    Returns:
        A tuple containing
        - the replacement (Expression)
        - the matched kernel (MatchedKernel)
        If no match was found, the replacement is the input expression and
        instead of the matched kernel, None is returned.
    """
    for node, pos in expr.preorder_iter():
        for kernel, substitution in discrimination_net.match(node):
            matched_kernel = kernel.set_match(substitution, False)
            expr = matchpy.replace(expr, pos, matched_kernel.replacement)
            return expr, matched_kernel
    return (expr, None)


def select_optimal_match(matches):
    """Selects the best match according to the cost function of the kernels.

    Selects the best match from matches according to the cost function. The
    matching kernels have to be equivalent in the sense that the patterns
    (excluding constraints) have to be the same (modulo renaming of variables).
    Lists of matches that fulfill those requirements can be obtained from
    ManyToOneMatcher.match(expr).grouped() (all matches in one group are
    equivalent) or with list(DiscriminationNet.match(expr)) (since
    discrimination nets only work for syntactic patterns, if there are multiple
    matches, all matching patterns have the same structure).

    Furthermore, the constraints of those patterns have to be lists of
    PropertyConstraint objects.

    The best matching kernel is selected in two steps. In the first step,
    the partial ordering on PropertyTuple objects is used to exclude all matches
    that are not minimal in the sense that there is another kernel that requires
    operands to have more specific properties. The underlying assumption is that
    more specialized kernels are expected to be cheaper. After this step, more
    than one kernel might be left because some kernels are not comparable.

    In the next step, the cost function is used to select the cheapest kernel
    among the remaining kernels.

    Args:
        matches (iterable): Iterable of tuples (Kernel, Substitution).

    Returns:
        Kernel: The cheapest matching kernel. If no match was found, it is None.
        Substitution: The corresponding substitution for the kernel. If no match
            was found, it is None.
    """

    matches = list(matches)

    if len(matches) == 1:
        return matches[0]

    optimal_match = (None, None)
    min_cost = math.inf

    # TODO: possible "optimization" if cost function is very expensive:
    # Construct list of candidates first. If there is only one, don't compute
    # cost and return it directly.

    for match in matches:
        if not any(match2[0].property_tuple < match[0].property_tuple for match2 in matches):
            cost = match[0].cost(match[1])
            if cost < min_cost:
                optimal_match = match
                min_cost = cost

    return optimal_match