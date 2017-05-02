

from ..derivation.graph.constructive import DerivationGraph

from collections import namedtuple

import time
from .. import examples

from .. import config

from .. import temporaries
from ..derivation import partitioning

to_file_format = "{0} {1:.4f} {2}"
time_format = "{:.4f}"
print_format = "{0} {1:.4f} ({2})"


all_examples= [ examples.Example001,
                examples.Example002,
                examples.Example003,
                examples.Example004,
                # examples.Example005,
                examples.Example006,
                examples.Example007,
                examples.Example008,
                examples.Example009,
                examples.Example010,
                examples.Example011,
                examples.Example012,
                examples.Example013,
                examples.Example014,
                examples.Example015,
                examples.Example016,
                examples.Example017,
                examples.Example018,
                examples.Example019,
                examples.Example020,
                examples.Example021,
                examples.Example022,
                examples.Example023,
                examples.Example024,
                examples.Example025,
                examples.Example026,
                examples.Example027,
                examples.Example028,
                examples.Example029,
                examples.Example030,
                examples.Example031,
                examples.Example032,
                examples.Example033,
                examples.Example034,
                examples.Example035,
                examples.Example036,
                examples.Example037,
                examples.Example038,
                examples.Example039,
                examples.Example040,
                examples.Example041,
                examples.Example042,
                examples.Example043,
                examples.Example044,
                examples.Example045,
                examples.Example046,
                examples.Example047,
                examples.Example048,
                examples.Example049,
                examples.Example050,
                examples.Example051,
                examples.Example052,
                examples.Example053,
                examples.Example054,
                examples.Example055,
                examples.Example056,
                examples.Example057,
                examples.Example058,
                examples.Example059,
                examples.Example060,
                examples.Example061,
                examples.Example062,
                examples.Example063,
                examples.Example064,
                examples.Example065,
                examples.Example066,
                examples.Example067,
                examples.Example068,
                examples.Example069,
                examples.Example070,
                examples.Example071,
                examples.Example072,
                examples.Example073,
                examples.Example074,
                examples.Example075,
                examples.Example076,
                examples.Example077,
                examples.Example078,
                examples.Example079,
                examples.Example080,
                examples.Example081,
                examples.Example082,
                # examples.Example083, # code generation: operand not in memory
                examples.Example084,
                examples.Example085,
                examples.Example086,
                examples.Example087,
                examples.Example088,
                examples.Example089,
                # examples.Example090, # does not work becaue of BlockedExpression
                # examples.Example091, # does not work becaue of BlockedExpression
                examples.Example092,
                examples.Example093,
                examples.Example094,
                examples.Example095,
                examples.Example096,
                examples.Example097,
                examples.Example098,
                examples.Example099,
                examples.Example100,
                examples.Example101,
                examples.Example102,
                examples.Example103,
                examples.Example104,
                examples.Example105,
                examples.Example106,
                examples.Example107,
                examples.Example108,
                examples.Example109,
                examples.Example110,
                examples.Example111,
                examples.Example112,]

def run(save_output=False):

    output = []

    output_file_path = "".join([config.language.name, "/example_trace.txt"])
    if save_output:
        output_file = open(output_file_path, "w")
    else:
        output_file = open(output_file_path)

    for example in all_examples:
        name = example.__name__
        graph = DerivationGraph(example().eqns)
        t_before = time.process_time()
        trace = graph.derivation(10, 6, False)
        execution_time = time.process_time() - t_before
        # print(example().eqns.to_julia_expression())
        if save_output:
            print(name, time_format.format(execution_time))
            print(to_file_format.format(name, execution_time, ":".join(str(t) for t in trace)), file=output_file)
            graph.to_dot_file("".join([config.language.name, "/graphs/", example.__name__, ".gv"]))

        else:
            graph.algorithms_to_files(pseudocode=True)
            _, old_execution_time, old_trace = next(output_file).split(" ")
            old_trace = old_trace.rstrip() # removing newline character
            print_string = print_format.format(name, execution_time, old_execution_time)
            # print(":".join(str(t) for t in trace))
            if old_trace != ":".join(str(t) for t in trace):
                print(print_string, "changed")
            else:
                print(print_string)
        temporaries.clear()
        partitioning.clear()

    output_file.close()