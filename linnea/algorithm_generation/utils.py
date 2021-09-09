import math


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