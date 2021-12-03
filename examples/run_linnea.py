



if __name__ == "__main__":

    import linnea.config

    linnea.config.set_output_code_path(".")
    linnea.config.init()

    from linnea.data_access.database import connect
    from linnea.algorithm_generation.graph.search_graph import SearchGraph

    from input1 import equations

    # import linnea.examples.examples
    # equations = linnea.examples.examples.Example001().eqns

    CREATE_ALGORITHMS_TABLE = "CREATE TABLE IF NOT EXISTS algorithms (id INTEGER PRIMARY KEY, equation TEXT, algorithm TEXT);"
    INSERT_ALGORITHM = "INSERT INTO algorithms (equation, algorithm) VALUES (?, ?);"


    # connect = sqlite3.connect("data.db")
    # with connect:
    #     connect.execute(CREATE_ALGORITHMS_TABLE)


    graph = SearchGraph(equations)

    graph.generate(time_limit=60,
                   merging=True,
                   dead_ends=True,
                   pruning_factor=1.0)

    graph.write_output(code=True,
                       generation_steps=False,
                       output_name="tmp",
                       experiment_code=False,
                       algorithms_limit=100,
                       graph=False)

    output = [(algorithm.code_as_function(), algorithm.generation_steps_website())  for algorithm in graph.k_best_algorithms(1)]
    
    # with connect:
    #     connect.execute(INSERT_ALGORITHM,  (output, equations))

    # connect.close() 