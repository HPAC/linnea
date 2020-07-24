from .. import temporaries

def add_expression(expr, properties):

    tmp = temporaries.create_tmp(expr)
    for prop in properties:
        tmp.set_property(prop)
