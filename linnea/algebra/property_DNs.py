
import matchpy

from .properties import Property as properties

from . import expression as ae

from .. import utils as lu

SPD_DN = None

def _init():
    """
    Due to import dependencies, these DNs can not be constructed during the
    initial imports. Instead, it is necessary to call this function before
    Linnea can be used.
    """

    SY1 = matchpy.Wildcard.symbol("SY1")
    WD1 = matchpy.Wildcard.dot("WD1")

    spd_constraint1 = lambda SY1: SY1.has_property(properties.FULL_RANK) and \
                                (SY1.has_property(properties.COLUMN_PANEL) or \
                                 SY1.has_property(properties.SQUARE))

    spd_constraint2 = lambda SY1, WD1: SY1.has_property(properties.VECTOR) or \
                            (SY1.has_property(properties.FULL_RANK) and \
                              (SY1.has_property(properties.COLUMN_PANEL) or \
                               SY1.has_property(properties.SQUARE)
                              )
                            ) and \
                            WD1.has_property(properties.SPD)

    spd_constraint3 = lambda SY1: SY1.has_property(properties.FULL_RANK) and \
                                (SY1.has_property(properties.ROW_PANEL) or \
                                 SY1.has_property(properties.SQUARE))

    spd_constraint4 = lambda SY1, WD1: SY1.has_property(properties.VECTOR) or \
                            (SY1.has_property(properties.FULL_RANK) and \
                              (SY1.has_property(properties.ROW_PANEL) or \
                               SY1.has_property(properties.SQUARE)
                              )
                            ) and \
                            WD1.has_property(properties.SPD)

    spd_pattern1 = matchpy.Pattern(ae.Times(ae.Transpose(SY1), SY1), matchpy.CustomConstraint(spd_constraint1))
    # spd_pattern5 = matchpy.Pattern(ae.Times(ae.InverseTranspose(SY1), ae.Inverse(SY1)), matchpy.CustomConstraint(spd_constraint1))
    spd_pattern2 = matchpy.Pattern(ae.Times(ae.Transpose(SY1), WD1, SY1), matchpy.CustomConstraint(spd_constraint2))
    spd_pattern3 = matchpy.Pattern(ae.Times(SY1, ae.Transpose(SY1)), matchpy.CustomConstraint(spd_constraint3))
    spd_pattern4 = matchpy.Pattern(ae.Times(SY1, WD1, ae.Transpose(SY1)), matchpy.CustomConstraint(spd_constraint4))

    global SPD_DN
    SPD_DN = matchpy.DiscriminationNet()
    SPD_DN.add(spd_pattern1)
    # SPD_DN.add(spd_pattern5)
    SPD_DN.add(spd_pattern2)
    SPD_DN.add(spd_pattern3)
    SPD_DN.add(spd_pattern4)


    spsd_constraint = lu.PropertyConstraint(WD1.variable_name, {properties.SPSD})

    spsd_pattern1 = matchpy.Pattern(ae.Times(ae.Transpose(SY1), SY1))
    spsd_pattern2 = matchpy.Pattern(ae.Times(ae.Transpose(SY1), WD1, SY1), spsd_constraint)
    spsd_pattern3 = matchpy.Pattern(ae.Times(SY1, ae.Transpose(SY1)))
    spsd_pattern4 = matchpy.Pattern(ae.Times(SY1, WD1, ae.Transpose(SY1)), spsd_constraint)

    global SPSD_DN
    SPSD_DN = matchpy.DiscriminationNet()
    SPSD_DN.add(spsd_pattern1)
    SPSD_DN.add(spsd_pattern2)
    SPSD_DN.add(spsd_pattern3)
    SPSD_DN.add(spsd_pattern4)