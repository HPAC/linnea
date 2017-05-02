
import matchpy

from .properties import Property as properties

from . import expression as ae

SPD_DN = None

def _init():
    """
    Due to import dependencies, these DNs can not be constructed during the
    initial imports. Instead, it is necessary to call this function before
    Linnea can be used.
    """

    SY1 = matchpy.Wildcard.symbol("SY1")
    WD1 = matchpy.Wildcard.dot("WD1")

    constraint1 = lambda SY1: SY1.has_property(properties.FULL_RANK) and \
                                (SY1.has_property(properties.COLUMN_PANEL) or \
                                 SY1.has_property(properties.SQUARE))

    constraint2 = lambda SY1, WD1: SY1.has_property(properties.VECTOR) or \
                            (SY1.has_property(properties.FULL_RANK) and \
                              (SY1.has_property(properties.COLUMN_PANEL) or \
                               SY1.has_property(properties.SQUARE)
                              )
                            ) and \
                            WD1.has_property(properties.SPD)

    constraint3 = lambda SY1: SY1.has_property(properties.FULL_RANK) and \
                                (SY1.has_property(properties.ROW_PANEL) or \
                                 SY1.has_property(properties.SQUARE))

    constraint4 = lambda SY1, WD1: SY1.has_property(properties.VECTOR) or \
                            (SY1.has_property(properties.FULL_RANK) and \
                              (SY1.has_property(properties.ROW_PANEL) or \
                               SY1.has_property(properties.SQUARE)
                              )
                            ) and \
                            WD1.has_property(properties.SPD)

    spd_pattern1 = matchpy.Pattern(ae.Times(ae.Transpose(SY1), SY1), matchpy.CustomConstraint(constraint1))
    spd_pattern2 = matchpy.Pattern(ae.Times(ae.Transpose(SY1), WD1, SY1), matchpy.CustomConstraint(constraint2))
    spd_pattern3 = matchpy.Pattern(ae.Times(SY1, ae.Transpose(SY1)), matchpy.CustomConstraint(constraint3))
    spd_pattern4 = matchpy.Pattern(ae.Times(SY1, WD1, ae.Transpose(SY1)), matchpy.CustomConstraint(constraint4))

    global SPD_DN
    SPD_DN = matchpy.DiscriminationNet()
    SPD_DN.add(spd_pattern1)
    SPD_DN.add(spd_pattern2)
    SPD_DN.add(spd_pattern3)
    SPD_DN.add(spd_pattern4)