
if __name__ == "__main__":

    import linnea.config

    linnea.config.set_output_code_path(".")
    linnea.config.init()

    from linnea.derivation.graph.derivation import DerivationGraph

    from input1 import equations

    # import linnea.examples.examples
    # equations = linnea.examples.examples.Example001().eqns

    graph = DerivationGraph(equations)
    graph.derivation(time_limit=60,
                     merging=True,
                     dead_ends=True,
                     pruning_factor=1.0)

    graph.write_output(code=True,
                       derivation=False,
                       output_name="tmp",
                       experiment_code=False,
                       algorithms_limit=100,
                       graph=False)