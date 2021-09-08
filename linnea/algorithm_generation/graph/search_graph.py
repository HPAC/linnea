from ...algebra import expression as ae

from ...utils import roundrobin

from .. import special_properties

from .utils import GenerationStep, generate_representations

from .base import equations

import itertools
import functools

class SearchGraph(equations.EquationsGraph):

    def generate(self, time_limit=60, merging=True, dead_ends=True, pruning_factor=1.0):
        
        # TODO add argument for stopping as soon as first solution is found
        # or use time_limit == 0?

        self.root.equations.infer_missing_properties()
        self.root.equations.infer_lhs_properties()
        self.root.equations.check_validity()
        self.root.equations = self.root.equations.to_normalform()

        trace_data, terminal_nodes = self.best_first_search(time_limit=time_limit, merging=merging, dead_ends=dead_ends, pruning_factor=pruning_factor)
                
        # file_name = os.path.join(config.output_code_path, config.output_name, config.language.name, "dfs_trace.csv")
        # with open(file_name, "wt", encoding='utf-8') as output_file:
        #     output_file.write("time, cost\n")
        #     for t, cost in trace_data:
        #         output_file.write("".join([str(t), ", ", str(cost), "\n"]))

        self.print_line()
        self.print_result("Number of nodes:", len(self.nodes))
        self.print_result("Solution nodes:", len(terminal_nodes))

        data = self.root.equations.get_data()
        self.print_result("Data:", "{:.3g}".format(data))
        if terminal_nodes:
            _, cost = self.shortest_path()
            self.print_result("Best solution:", "{:.3g}".format(cost))
            self.print_result("Intensity:", "{:.3g}".format(cost/data))        

        # from ... import temporaries
        # for tmp, expr in temporaries._equivalent_expressions.items():
        #     print(tmp, "<->", expr)
        #     for _expr, _tmp in temporaries._table_of_temporaries.items():
        #         if tmp == _tmp.name and expr != _expr:
        #             print("<-", _expr)

        return trace_data


    def successor_generator(self, node):

        equation, eqn_idx = node.equations.process_next()
        if not equation:
            eqn_idx = None

        gen_var_single1, gen_var_single2, gen_var_single3 = itertools.tee(generate_representations(node.equations, eqn_idx), 3)
        gen_var1, gen_var2 = itertools.tee(generate_representations(node.equations), 2)

        kernels_constructive = map(functools.partial(self.GS_kernels_constructive, node), gen_var_single1)
        kernels = map(functools.partial(self.GS_kernels, node), gen_var_single2)
        CSE_replacement = map(functools.partial(self.GS_CSE_replacement, node), gen_var1)
        tricks = map(functools.partial(self.GS_tricks, node), gen_var_single3)
        factorizations = map(functools.partial(self.GS_factorizations, node), gen_var2)

        if GenerationStep.factorizations in node.applied_generation_steps:
            funs = [kernels_constructive, CSE_replacement, kernels, tricks]
        elif GenerationStep.kernels in node.applied_generation_steps:
            funs = [kernels_constructive, CSE_replacement, factorizations, kernels, tricks]
        elif GenerationStep.CSE in node.applied_generation_steps:
            # there are cases where doing CSE replacement twice results in better solutions
            funs = [kernels_constructive, CSE_replacement, factorizations, kernels, tricks]
        elif GenerationStep.tricks in node.applied_generation_steps:
            funs = [kernels_constructive, CSE_replacement, factorizations, kernels]
        else:
            funs = [CSE_replacement, kernels_constructive, factorizations, kernels, tricks]

        yield from roundrobin(*roundrobin(*funs))


    # def successor_generator(self, node):

    #     if GenerationStep.factorizations in node.applied_generation_steps:
    #         funs = [self.GS_kernels_constructive, self.GS_kernels, self.GS_CSE_replacement, self.GS_tricks]
    #     elif GenerationStep.kernels in node.applied_generation_steps:
    #         funs = [self.GS_kernels_constructive, self.GS_kernels, self.GS_CSE_replacement, self.GS_factorizations, self.GS_tricks]
    #     elif GenerationStep.CSE in node.applied_generation_steps:
    #         # there are cases where doing CSE replacement twice results in better solutions
    #         funs = [self.GS_kernels_constructive, self.GS_kernels, self.GS_CSE_replacement, self.GS_factorizations, self.GS_tricks]
    #     elif GenerationStep.tricks in node.applied_generation_steps:
    #         funs = [self.GS_kernels_constructive, self.GS_kernels, self.GS_CSE_replacement, self.GS_factorizations]
    #     else:
    #         funs = [self.GS_CSE_replacement, self.GS_kernels_constructive, self.GS_kernels, self.GS_factorizations, self.GS_tricks]

    #     generators = [fun(node, eqns) for eqns, fun in itertools.product(generate_representations(node.equations), funs)]
    #     yield from self.roundrobin(*generators)
