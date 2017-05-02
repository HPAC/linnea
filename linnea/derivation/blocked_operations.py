
# Used to store blocked/forbidden products. Usually, forbidden products
# are factors obtained from the same factorization.

_blocked = set()

def is_blocked(op):
    # expects a string of the expression that represents the operations
    # TODO do I really want that?
    global _blocked
    return op in _blocked

def set_blocked(ops):
    global _blocked
    _blocked |= set(ops)


