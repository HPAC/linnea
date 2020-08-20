
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
        utils.CodeTemplate(textwrap.dedent(
            """\
            $output = zeros($type, $m, $n)
            for i = 1:min($m, $n);
                $output[i, i] = $input[i];
            end;
            """
            ))
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
        utils.CodeTemplate("$output = map(x -> x[2], sort!(findall(!iszero, $input), by = x -> x[1]))\n") # starting from Julia 1.1, this could also be done with $output = invperm(map(x -> findall(!iszero, x)[1], eachcol($input)))
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
            end;
            """
            ))
        ),
    sf.StorageFormatConversion(
        sf.StorageFormat.ipiv,
        sf.StorageFormat.full,
        utils.CodeTemplate(textwrap.dedent(
            """\
            tmp = [1:length($input);]
            @inbounds for i in 1:length($input)
                tmp[i], tmp[$input[i]] = tmp[$input[i]], tmp[i];
            end;
            $output = Array{$type}(I, $n, $n)[tmp,:]
            """
            ))
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
        utils.CodeTemplate(textwrap.dedent(
            """\
            triu!($op, 1)
            for i = 1:min(size($op)...);
                $op[i,i] = one($type);
            end;
            """
            ))
        ),
    sf.StorageFormatConversion(
        sf.StorageFormat.lower_triangular_udiag,
        sf.StorageFormat.full,
        utils.CodeTemplate(textwrap.dedent(
            """\
            tril!($op, -1)
            for i = 1:min(size($op)...);
                $op[i,i] = one($type);
            end;
            """
            ))
        ),
    sf.StorageFormatConversion(
        sf.StorageFormat.symmetric_upper_triangular,
        sf.StorageFormat.full,
        utils.CodeTemplate(textwrap.dedent(
            """\
            for i = 1:$n-1;
                view($op, i+1:$n, i)[:] = view($op, i, i+1:$n);
            end;
            """
            ))
        ),
    sf.StorageFormatConversion(
        sf.StorageFormat.symmetric_lower_triangular,
        sf.StorageFormat.full,
        utils.CodeTemplate(textwrap.dedent(
            """\
            for i = 1:$n-1;
                view($op, i, i+1:$n)[:] = view($op, i+1:$n, i);
            end;
            """
            ))
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