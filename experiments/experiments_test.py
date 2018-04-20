
if __name__ == "__main__":

    import linnea.config

    linnea.config.init()
    linnea.config.set_output_path("~/linnea/output/")

    from linnea.derivation.graph.constructive import DerivationGraph
    # from linnea.derivation.graph.exhaustive import DerivationGraph
    # from linnea.derivation.graph.matrix_chain_derivation import DerivationGraph

    import linnea.examples.examples
    import linnea.examples.lamp_paper.examples

    example = linnea.examples.examples.Example063()

    # print(example.eqns)
    graph = DerivationGraph(example.eqns)
    trace = graph.derivation(
                        solution_nodes_limit=100,
                        iteration_limit=15,
                        merging=True,
                        dead_ends=True)
    print(":".join(str(t) for t in trace))

    # import linnea.temporaries
    # print("\n".join(["{}: {}".format(k, v) for k, v in linnea.temporaries._table_of_temporaries.items()]))
    # print("\n".join(["{}: {}".format(k, v) for k, v in linnea.temporaries._equivalent_expressions.items()]))

    graph.write_output(code=True,
                       pseudocode=True,
                       output_name="tmp",
                       operand_generator=True,
                       algorithms_limit=100,
                       graph=True)



