
if __name__ == "__main__":

    import linnea.config

    linnea.config.set_output_path(".")
    linnea.config.init()

    from linnea.derivation.graph.constructive import DerivationGraph
    # from linnea.derivation.graph.exhaustive import DerivationGraph

    from input1 import equations

    # import linnea.examples.examples
    # equations = linnea.examples.examples.Example001().eqns

    graph = DerivationGraph(equations)
    graph.derivation(solution_nodes_limit=100,
                     iteration_limit=100,
                     merging=True)

    graph.write_output(code=True,
                       pseudocode=False,
                       output_name="tmp",
                       operand_generator=False,
                       algorithms_limit=100,
                       graph=False)