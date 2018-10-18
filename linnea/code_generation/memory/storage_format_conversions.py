
from  . import storage_format as sf

from ... import utils

import itertools
import operator
import textwrap

out_of_place_conversions = [
    sf.StorageFormatConversion(
        sf.StorageFormat.full,
        sf.StorageFormat.diagonal_vector,
        utils.CodeTemplate("$output = diag($input)\n")
        ),
    sf.StorageFormatConversion(
        sf.StorageFormat.diagonal_vector,
        sf.StorageFormat.full,
        # utils.CodeTemplate("$output = diagm($input)\n")
        utils.CodeTemplate("$output = zeros($type, $m, $n)\nfor i = 1:min($m, $n); $output[i, i] = $input[i]; end\n")
        ),
    sf.StorageFormatConversion(
        sf.StorageFormat.symmetric_lower_triangular,
        sf.StorageFormat.diagonal_vector,
        utils.CodeTemplate("$output = diag($input)\n")
        ),
    sf.StorageFormatConversion(
        sf.StorageFormat.symmetric_upper_triangular,
        sf.StorageFormat.diagonal_vector,
        utils.CodeTemplate("$output = diag($input)\n")
        ),
    sf.StorageFormatConversion(
        sf.StorageFormat.full,
        sf.StorageFormat.permutation_vector,
        utils.CodeTemplate("$output = invperm(findn($input)[1])\n") # this could also be done with convert(Array{Int64, 1}, P*range(1, size(P, 1)))
        ),
    sf.StorageFormatConversion(
        sf.StorageFormat.permutation_vector,
        sf.StorageFormat.full,
        utils.CodeTemplate("$output = Array{$type}(I, $n, $n)[$input,:]\n")
        ),
    sf.StorageFormatConversion(
        sf.StorageFormat.ipiv,
        sf.StorageFormat.permutation_vector,
        utils.CodeTemplate(textwrap.dedent(
            """\
            $output = [1:length($input);]
            @inbounds for i in 1:length($input)
                $output[i], $output[$input[i]] = $output[$input[i]], $output[i];
            end
            """
            ))
        ),
    sf.StorageFormatConversion(
        sf.StorageFormat.cholfact_L,
        sf.StorageFormat.full,
        utils.CodeTemplate("$output = convert(Array{$type, 2}, $input[:L])\n")
        ),
    sf.StorageFormatConversion(
        sf.StorageFormat.LUfact_L,
        sf.StorageFormat.full,
        utils.CodeTemplate("$output = $input[:L]\n")
        ),
    sf.StorageFormatConversion(
        sf.StorageFormat.LUfact_U,
        sf.StorageFormat.full,
        utils.CodeTemplate("$output = $input[:U]\n")
        ),
    sf.StorageFormatConversion(
        sf.StorageFormat.LUfact_P,
        sf.StorageFormat.full,
        utils.CodeTemplate("$output = $input[:P]\n")
        ),
    sf.StorageFormatConversion(
        sf.StorageFormat.LUfact_P,
        sf.StorageFormat.permutation_vector,
        utils.CodeTemplate("$output = $input[:p]\n")
        ),
    sf.StorageFormatConversion(
        sf.StorageFormat.QRfact_Q,
        sf.StorageFormat.full,
        utils.CodeTemplate("$output = Array($input.Q)\n")
        ),
    sf.StorageFormatConversion(
        sf.StorageFormat.QRfact_R,
        sf.StorageFormat.full,
        utils.CodeTemplate("$output = $input.R\n")
        ),
    sf.StorageFormatConversion(
        sf.StorageFormat.eigfact_Z,
        sf.StorageFormat.full,
        utils.CodeTemplate("$output = $input[:vectors]\n")
        ),
    sf.StorageFormatConversion(
        sf.StorageFormat.eigfact_W,
        sf.StorageFormat.diagonal_vector,
        utils.CodeTemplate("$output = $input[:values]\n")
        ),
    sf.StorageFormatConversion(
        sf.StorageFormat.eigfact_W,
        sf.StorageFormat.full,
        utils.CodeTemplate("$output = diagm($input[:values])\n")
        ),
    sf.StorageFormatConversion(
        sf.StorageFormat.svdfact_U,
        sf.StorageFormat.full,
        utils.CodeTemplate("$output = $input.U\n")
        ),
    sf.StorageFormatConversion(
        sf.StorageFormat.svdfact_S,
        sf.StorageFormat.diagonal_vector,
        utils.CodeTemplate("$output = $input.S\n")
        ),
    sf.StorageFormatConversion(
        sf.StorageFormat.svdfact_S,
        sf.StorageFormat.full,
        utils.CodeTemplate("$output = zeros($type, $m, $n)\nd = $input.S\nfor i = 1:min($m, $n); $output[i, i] = d[i]; end\n")
        ),
    sf.StorageFormatConversion(
        sf.StorageFormat.svdfact_V,
        sf.StorageFormat.full,
        utils.CodeTemplate("$output = $input.Vt\n")
        ),
]

in_place_conversions = [
    sf.StorageFormatConversion(
        sf.StorageFormat.upper_triangular,
        sf.StorageFormat.full,
        utils.CodeTemplate("triu!($op)\n")
        ),
    sf.StorageFormatConversion(
        sf.StorageFormat.lower_triangular,
        sf.StorageFormat.full,
        utils.CodeTemplate("tril!($op)\n")
        ),
    sf.StorageFormatConversion(
        sf.StorageFormat.upper_triangular_udiag,
        sf.StorageFormat.full,
        utils.CodeTemplate("triu!($op, 1)\nfor i = 1:min(size($op)...); $op[i,i] = one($type); end\n")
        ),
    sf.StorageFormatConversion(
        sf.StorageFormat.lower_triangular_udiag,
        sf.StorageFormat.full,
        utils.CodeTemplate("tril!($op, -1)\nfor i = 1:min(size($op)...); $op[i,i] = one($type); end\n")
        ),
    sf.StorageFormatConversion(
        sf.StorageFormat.symmetric_upper_triangular,
        sf.StorageFormat.full,
        utils.CodeTemplate("for i = 1:$n-1; $op[i+1:$n,i] = $op[i,i+1:$n]; end;\n")
        ),
    sf.StorageFormatConversion(
        sf.StorageFormat.symmetric_lower_triangular,
        sf.StorageFormat.full,
        utils.CodeTemplate("for i = 1:$n-1; $op[i,i+1:$n] = $op[i+1:$n,i]; end;\n")
        ),
    ]

def _construct_conversion_dict(conversions):

    conversions.sort(key=operator.attrgetter("source_format.value"))

    conversion_dict = dict()
    for key, group in itertools.groupby(conversions, operator.attrgetter("source_format")):
        # print(key, list(group))
        conversion_dict[key] = list(group)

    return conversion_dict

in_place_conversion_dict = _construct_conversion_dict(in_place_conversions)
out_of_place_conversion_dict = _construct_conversion_dict(out_of_place_conversions)


if __name__ == "__main__":

    d = _construct_conversion_dict(out_of_place_conversions)