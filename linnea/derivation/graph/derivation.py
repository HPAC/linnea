from ...algebra import expression as ae

from ...utils import is_inverse

from .. import special_properties

from .utils import DS_step, generate_variants, find_operands_to_factor, \
                   find_explicit_symbol_inverse

from . import base

import itertools

class DerivationGraph(base.derivation.DerivationGraphBase):

    def derivation(self, time_limit=60, merging=True, dead_ends=True):
        
        # TODO add argument for stopping as soon as first solution is found
        # or use time_limit == 0?

        self.root.equations.check_validity()
        self.root.equations = self.root.equations.to_normalform()
        self.root.equations.infer_lhs_properties()

        self.init_temporaries(self.root.equations)
        self.root.equations.check_consistency()

        trace_data, terminal_nodes = self.best_first_search(time_limit=time_limit, merging=merging, dead_ends=dead_ends)
                
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

        # from .... import temporaries
        # print("\n".join(["{}: {}".format(k, v) for k, v in temporaries._equivalent_expressions.items()]))
        # print("######")
        # print("\n".join(["{}: {}".format(k, v) for k, v in temporaries._table_of_temporaries.items()]))

        return trace_data


    def successor_generator(self, node):

        _, eqn_idx = node.equations.process_next()
        
        gen_var_single1, gen_var_single2, gen_var_single3 = itertools.tee(generate_variants(node.equations, eqn_idx), 3)
        gen_var1, gen_var2 = itertools.tee(generate_variants(node.equations), 2)

        import functools
        kernels_constructive = map(functools.partial(self.DFS_kernels_constructive, node), gen_var_single1)
        kernels = map(functools.partial(self.DFS_kernels, node), gen_var_single2)
        CSE_replacement = map(functools.partial(self.DFS_CSE_replacement, node), gen_var1)
        tricks = map(functools.partial(self.DFS_tricks, node), gen_var_single3)
        factorizations = map(functools.partial(self.DFS_factorizations, node), gen_var2)

        if DS_step.factorizations in node.applied_DS_steps:
            funs = [kernels_constructive, kernels, CSE_replacement, tricks]
        elif DS_step.kernels in node.applied_DS_steps:
            funs = [kernels_constructive, kernels, CSE_replacement, factorizations, tricks]
        elif DS_step.CSE in node.applied_DS_steps:
            # there are cases where doing CSE replacement twice results in better solutions
            funs = [kernels_constructive, kernels, CSE_replacement, factorizations, tricks]
        elif DS_step.tricks in node.applied_DS_steps:
            funs = [kernels_constructive, kernels, CSE_replacement, factorizations]
        else:
            funs = [CSE_replacement, kernels_constructive, kernels, factorizations, tricks]

        yield from self.roundrobin(*self.roundrobin(*funs))


    def successor_generator_(self, node):

        if DS_step.factorizations in node.applied_DS_steps:
            funs = [self.DFS_kernels_constructive, self.DFS_kernels, self.DFS_CSE_replacement, self.DFS_tricks]
        elif DS_step.kernels in node.applied_DS_steps:
            funs = [self.DFS_kernels_constructive, self.DFS_kernels, self.DFS_CSE_replacement, self.DFS_factorizations, self.DFS_tricks]
        elif DS_step.CSE in node.applied_DS_steps:
            # there are cases where doing CSE replacement twice results in better solutions
            funs = [self.DFS_kernels_constructive, self.DFS_kernels, self.DFS_CSE_replacement, self.DFS_factorizations, self.DFS_tricks]
        elif DS_step.tricks in node.applied_DS_steps:
            funs = [self.DFS_kernels_constructive, self.DFS_kernels, self.DFS_CSE_replacement, self.DFS_factorizations]
        else:
            funs = [self.DFS_CSE_replacement, self.DFS_kernels_constructive, self.DFS_kernels, self.DFS_factorizations, self.DFS_tricks]

        generators = [fun(node, eqns) for eqns, fun in itertools.product(generate_variants(node.equations), funs)]
        yield from self.roundrobin(*generators)


    def init_temporaries(self, equations):

        seen_before = set()
        operands_to_factor = find_operands_to_factor(equations)
        for equations_var in generate_variants(equations):
            for equation in equations_var:
                for inv_expr, pos in find_explicit_symbol_inverse(equation.rhs):
                    if inv_expr not in seen_before:
                        special_properties.add_expression(inv_expr, [])
                        seen_before.add(inv_expr)
                for expr, pos in equation.rhs.preorder_iter():
                    # if isinstance(expr, Times) and any(is_inverse(operand) and (not isinstance(operand.operand, Symbol) or operand.operand in operands_to_factor) for operand in expr.operands):
                    if (isinstance(expr, ae.Times) and any(is_inverse(operand) and operand.operand in operands_to_factor for operand in expr.operands)) or (is_inverse(expr) and not isinstance(expr.operand, ae.Symbol)):
                        # TODO is it dangerous to add inv(expr) here? It becomes explicit inversion, even if it originally wasn't.
                        if expr not in seen_before:
                            special_properties.add_expression(expr, [])
                            seen_before.add(expr)
