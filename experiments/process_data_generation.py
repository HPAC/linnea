import os.path
import pandas as pd

def read_results():
    # what about missing data?

    data = []
    for jobindex in range(1, 32):
        file_name = "linnea_generation{}.csv".format(jobindex)
        if os.path.isfile(file_name):
            data.append(pd.read_csv(file_name, index_col=[0, 1, 2]))
        else:
            print("Missing file", file_name)

    return pd.concat(data)

    

if __name__ == '__main__':

    results = read_results()
    results.to_csv("linnea_generation.csv")
    # print(results)

    new_columns = ["constructive_merging", "constructive_no_merging", "exhaustive_merging", "exhaustive_no_merging"]

    # we have to do something about missing data here

    nodes_count = results.drop(columns=["mean", "std", "min", "max", "solution"])
    nodes_count = nodes_count.unstack([1, 2]) # use fill_value=... for missing data
    nodes_count.columns = new_columns
    # nodes_count.index.name = None
    nodes_count.to_csv("generation_nodes_count.csv")
    # print(nodes_count)

    generation_time = results.drop(columns=["std", "min", "max", "nodes", "solution"])
    generation_time = generation_time.unstack([1, 2]) # use fill_value=... for missing data
    generation_time.columns = new_columns
    # generation_time.index.name = None
    generation_time.to_csv("generation_time.csv")
    # print(generation_time)