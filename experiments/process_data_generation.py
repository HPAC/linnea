import os.path
import argparse
import pandas as pd

def read_results(number_of_experiments, dirs):
    # what about missing data?

    data = []
    for dir in dirs:
        for jobindex in range(1, number_of_experiments+1):
            file_name = "{}/linnea_generation{}.csv".format(dir, jobindex)
            if os.path.isfile(file_name):
                data.append(pd.read_csv(file_name, index_col=[0, 1, 2]))
            else:
                print("Missing file", file_name)

    return pd.concat(data)

    

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog="process_data")
    parser.add_argument("experiment", choices=["lamp_example", "random", "pldi"])
    args = parser.parse_args()

    if args.experiment == "lamp_example":
        number_of_experiments = 31
    elif args.experiment == "random":
        number_of_experiments = 100
    elif args.experiment == "pldi":
        number_of_experiments = 25

    dirs = []
    new_columns = []
    for dir, col in zip(["c_m", "c_nm", "e_m", "e_nm"], ["constructive_merging", "constructive_no_merging", "exhaustive_merging", "exhaustive_no_merging"]):
        if os.path.exists(dir):
            dirs.append(dir)
            new_columns.append(col)
        else:
            print("No directory", dir)

    results = read_results(number_of_experiments, dirs)
    results.to_csv("linnea_generation.csv", na_rep="NaN")

    # we have to do something about missing data here

    nodes_count = results.drop(columns=["mean", "std", "min", "max", "solution"])
    nodes_count = nodes_count.unstack([1, 2]) # use fill_value=... for missing data
    nodes_count.columns = new_columns
    # nodes_count.index.name = None
    nodes_count.to_csv("generation_nodes_count.csv", na_rep="NaN")
    # print(nodes_count)

    generation_time = results.drop(columns=["std", "min", "max", "nodes", "solution"])
    generation_time = generation_time.unstack([1, 2]) # use fill_value=... for missing data
    generation_time.columns = new_columns
    # generation_time.index.name = None
    generation_time.to_csv("generation_time.csv", na_rep="NaN")
    # print(generation_time)

    generation_time_renamed = generation_time.rename(mapper=lambda name: name + "_time", axis='columns')
    nodes_count_renamed = nodes_count.rename(mapper=lambda name: name + "_nodes", axis='columns')

    # all_names = ["constructive_merging_time", "constructive_no_merging_time", "exhaustive_merging_time", "exhaustive_no_merging_time", "constructive_merging_nodes", "constructive_no_merging_nodes", "exhaustive_merging_nodes", "exhaustive_no_merging_nodes"]
    nodes_and_time = pd.concat([generation_time_renamed, nodes_count_renamed], axis=1)
    # print(nodes_and_time)
    nodes_and_time.to_csv("generation_all.csv", na_rep="NaN")