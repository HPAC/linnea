
import enum
import copy

from ... import utils

@utils.PartiallyOrderedEnum
class StorageFormat(enum.Enum):
    full = 0
    upper_triangular = 1
    lower_triangular = 2
    diagonal_vector = 3
    permutation_vector = 4
    upper_triangular_udiag = 5
    lower_triangular_udiag = 6
    as_overwritten = 7 # output only. IMPORTANT: Do not use this for SYRK.
    symmetric_triangular = 8 # requirement only
    symmetric_upper_triangular = 9
    symmetric_lower_triangular = 10
    factorization_obj = 11
    triangular_udiag_opt = 12 # requirement only
    cholfact_L = 13
    LUfact_L = 14
    LUfact_U = 15
    LUfact_P = 16
    QRfact_Q = 17
    QRfact_R = 18
    eigfact_Z = 19
    eigfact_W = 20
    svdfact_U = 21
    svdfact_S = 22
    svdfact_V = 23
    as_vector = 24
    """
    Output only. Actual format will be either symmetric_upper_triangular or
    symmetric_lower_triangular. If both are possible, symmetric_lower_triangular
    will be chosen.
    """
    symmetric_triangular_out = 25
    ipiv = 26

    __ordering__ = {
        (triangular_udiag_opt, upper_triangular_udiag),
        (upper_triangular_udiag, upper_triangular),
        (upper_triangular, full),
        (triangular_udiag_opt, lower_triangular_udiag),
        (lower_triangular_udiag, lower_triangular),
        (lower_triangular, full),
        (symmetric_upper_triangular, full),
        (symmetric_lower_triangular, full),
        (symmetric_triangular, symmetric_upper_triangular),
        (symmetric_triangular, symmetric_lower_triangular),
        (cholfact_L, factorization_obj),
        (LUfact_L, factorization_obj),
        (LUfact_U, factorization_obj),
        (LUfact_P, factorization_obj),
        (QRfact_Q, factorization_obj),
        (QRfact_R, factorization_obj),
        (eigfact_Z, factorization_obj),
        (eigfact_W, factorization_obj),
        (svdfact_U, factorization_obj),
        (svdfact_S, factorization_obj),
        (svdfact_V, factorization_obj),
        (diagonal_vector, as_vector),
        (permutation_vector, as_vector),
        (ipiv, as_vector),
        (symmetric_lower_triangular, symmetric_triangular_out),
        (symmetric_upper_triangular, symmetric_triangular_out),
        }

"""
minimal set:
full
upper_triangular
lower_triangular
as_vector
implicit_unit_diagonal

as_overwritten # output only
triangular # requirement only

"""

class MissingStorageFormatConversion(Exception):
    pass

class IncompatibleStorageFormats(Exception):
    pass

def select_storage_format_conversion(source_format, target_format, conversion_dict):
    """Selects a suitable storage format converions.

    Selects a conversion from source_format to a target format that is
    compatible with target_format from the conversions stored in
    conversion_dict.  For example, if target_format is
    StorageFormat.lower_triangular, the selected conversion may convert to
    StorageFormat.full.

    Args:
        source_format (StorageFormat): The original storage format.
        target_format (StorageFormet): The format that the actual resulting
            format has to be compatible with.
        conversion_dict (dict): A dictionary that maps a StorageFormat object
            to a list of StorageFormatConversion objects. The key is the
            source_format of those conversions.

    Returns:
        StorageFormatConversion: A suitable storage format conversion, if there
            is one.

    Raises:
        MissingStorageFormatConversion: If no suitable StorageFormatConversion
            was found.
    """
    # print(source_format, target_format)
    for conversion in conversion_dict[source_format]:
        if target_format <= conversion.target_format:
            return conversion
    raise MissingStorageFormatConversion("No conversion from {0} to {1}.".format(source_format, target_format))
    # print("MISSING", source_format, target_format)

class StorageFormatConversion():
    """docstring for StorageFormatConversion"""
    def __init__(self, source_format, target_format, code_template):
        self.source_format = source_format
        self.target_format = target_format
        self.code_template = code_template

if __name__ == "__main__":
    pass
