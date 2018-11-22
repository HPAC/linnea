import os.path
import pandas as pd

def read_results_generation(experiment, number_of_experiments):
    # what about missing data?

    dirs = []
    new_columns = []
    for dir, col in zip(["c_m", "c_nm", "e_m", "e_nm"], ["constructive_merging", "constructive_no_merging", "exhaustive_merging", "exhaustive_no_merging"]):
        path = os.path.join(experiment, "generation", dir)
        if os.path.exists(path):
            dirs.append(path)
            new_columns.append(col)
        else:
            print("No directory", path)

    data = []
    for dir in dirs:
        for jobindex in range(1, number_of_experiments+1):
            file_name = "{}/generation{:03}.csv".format(dir, jobindex)
            if os.path.isfile(file_name):
                data.append(pd.read_csv(file_name, index_col=[0, 1, 2]))
            else:
                print("Missing file", file_name)

    if data:
        return pd.concat(data)
    else:
        return None


def read_results_execution(experiment_name, number_of_experiments, usecols=range(0, 2)):

    """ TODO idea: Instead of usecols, use Enum an argument, which is maped to
    columns. This can then also be used in df.rename()

    What if we want more than one column?
    """

    # TODO move this into the loop below
    example_names = ["{}{:03}".format(experiment_name, i) for i in range(1, number_of_experiments+1)]

    language_dfs = []
    for language in ["cpp", "julia", "matlab"]:
        path = os.path.join(experiment_name, "execution", language)
        if os.path.exists(path):
            file_dfs = []
            for example in example_names:
                file_path = "{0}/{1}_results_{2}.txt".format(path, language, example)
                if os.path.isfile(file_path):
                    df = pd.read_csv(file_path, sep='\t', skipinitialspace=True, usecols=usecols, index_col=0)
                    df = df.transpose()
                    # df.rename(index={"Time": example}, inplace=True)
                    # This works only when we extrace not more than one column
                    df.rename(mapper=lambda _: example, inplace=True)
                    file_dfs.append(df)
                else:
                    print("Missing file", file_path)


            df = pd.concat(file_dfs, sort=True)
            language_dfs.append(df)
        else:
            print("No results for", language)

    return pd.concat(language_dfs, axis=1, sort=True)


def to_performance_profiles_data(time_data):

    # normalize by fastest implementation
    normalized_data = time_data.apply(lambda row: row/row.min(), axis=1)

    # sort columns
    sorted_data = pd.concat([normalized_data[col].sort_values().reset_index(drop=True) for col in normalized_data], axis=1, ignore_index=True)
    sorted_data.columns = time_data.columns
    
    # first row (for plotting)
    ones = pd.DataFrame({col:1.0 for col in sorted_data.columns}, index=[0])

    # last row (for plotting)
    max_value = sorted_data.iloc[-1,:].max()*1.1
    max_values = pd.DataFrame({col:max_value for col in sorted_data.columns}, index=[0])

    # add first and last row
    sorted_data = pd.concat([ones, sorted_data, max_values], axis=0, ignore_index=True)

    # construct y axis
    n = len(sorted_data.index)

    y_coordinates = pd.Series([val/(n-2) for val in range(0, n)])
    y_coordinates.iloc[-1] = 1

    # add y axis
    final = pd.concat([y_coordinates, sorted_data], axis=1)
    final.rename(columns={0: "y"}, inplace=True)

    return final


def select_best(d1, d2):
    if d1 < d2:
        return d1
    else:
        return d2


def performance_profiles_data_reduce(time_data):

    # print(time_data)

    armadillo = time_data.apply(lambda row: select_best(row["naive_armadillo"], row["recommended_armadillo"]), axis=1)
    matlab = time_data.apply(lambda row: select_best(row["naive_matlab"], row["recommended_matlab"]), axis=1)
    julia = time_data.apply(lambda row: select_best(row["naive_julia"], row["recommended_julia"]), axis=1)
    eigen = time_data.apply(lambda row: select_best(row["naive_eigen"], row["recommended_eigen"]), axis=1)
    # print(armadillo)

    res = pd.concat([time_data["algorithm0e"], armadillo, matlab, julia, eigen], axis=1)
    res.rename(columns={0: "armadillo", 1: "matlab", 2: "julia", 3: "eigen"}, inplace=True)
    # print(res)

    return res


def to_speedup_data(time_data, reference):
    speedup_data = time_data.apply(lambda row: row/row[reference], axis=1)
    speedup_data["mean"] = speedup_data.apply(lambda row: row.mean(), axis=1)
    return speedup_data


def to_mean_speedup(speedup_data, drop_reference=None):
    speedup_mean = pd.DataFrame(speedup_data.mean())
    if drop_reference:
        speedup_mean.drop(drop_reference, inplace=True)
    speedup_mean = speedup_mean.transpose()
    return speedup_mean


def compare_to_fastest(time_data, reference):
    """Identifies cases where other systems are faster than reference.

    Returns a DataFrame with example problem, algorithm and slowdown, sorted by
    the magintude of the slowdown.
    """
    diff = time_data.apply(lambda row: row[reference]/row, axis=1)
    diff = diff.stack()
    diff = diff.loc[diff > 1]
    diff.sort_values(ascending=False, inplace=True)
    return(diff)


def check_std_dev(experiment_name, number_of_experiments, threshold=0.1):
    execution_time = read_results(experiment_name, number_of_experiments)
    std_dev = read_results(experiment_name, number_of_experiments, usecols=[0, 2])
    rel_std_dev = std_dev/execution_time.values
    rel_std_dev = rel_std_dev.stack()
    rel_std_dev = rel_std_dev.loc[rel_std_dev > threshold]
    rel_std_dev.sort_values(ascending=False, inplace=True)
    return rel_std_dev


def time_accumulated(time_data):
    time_mean = pd.DataFrame(time_data.mean(), columns=["mean"])
    time_min = pd.DataFrame(time_data.min(), columns=["min"])
    time_max = pd.DataFrame(time_data.max(), columns=["max"])
    
    accumulated = pd.concat([time_mean, time_min, time_max], axis=1)
    return accumulated


def get_column_names(experiment):

    if experiment == "combined":
        experiments = ["random", "pldi"]
    else:    
        experiments = [experiment]

    # TODO this code also appears in "read_results_generation"    
    new_columns = []
    for exp in experiments:
        for dir, col in zip(["c_m", "c_nm", "e_m", "e_nm"], ["constructive_merging", "constructive_no_merging", "exhaustive_merging", "exhaustive_no_merging"]):
            if col in new_columns:
                continue
            path = os.path.join(exp, "generation", dir)
            if os.path.exists(path):
                new_columns.append(col)

    return new_columns


def process_data_execution(execution_time, experiment):

    execution_time.to_csv("{}_execution_time.csv".format(experiment), na_rep="NaN")
    execution_time.dropna().to_csv("{}_execution_time_clean.csv".format(experiment))

    performance_profiles_data = to_performance_profiles_data(performance_profiles_data_reduce(execution_time))
    performance_profiles_data.to_csv("{}_performance_profile.csv".format(experiment), na_rep="NaN")
    
    speedup_data = to_speedup_data(execution_time, speedup_reference)
    speedup_data.to_csv("{}_speedup.csv".format(experiment), na_rep="NaN")
    speedup_data.dropna().to_csv("{}_speedup_clean.csv".format(experiment))

    mean_speedup = to_mean_speedup(speedup_data)
    mean_speedup.to_csv("{}_mean_speedup.csv".format(experiment), na_rep="NaN")

    speedup_over_linnea = compare_to_fastest(execution_time, speedup_reference)
    speedup_over_linnea.to_csv("{}_speedup_over_linnea.csv".format(experiment), na_rep="NaN")

def process_data_generation(gen_results, experiment):

    gen_results.to_csv("{}_generation.csv".format(experiment), na_rep="NaN")

    new_columns = get_column_names(experiment)
    nodes_count = gen_results.drop(columns=["mean", "std", "min", "max", "solution"])
    nodes_count = nodes_count.unstack([1, 2]) # use fill_value=... for missing data
    nodes_count.columns = new_columns
    nodes_count.to_csv("{}_generation_nodes_count.csv".format(experiment), na_rep="NaN")
    nodes_count.dropna().to_csv("{}_generation_nodes_count_clean.csv".format(experiment))

    generation_time = gen_results.drop(columns=["std", "min", "max", "nodes", "solution"])
    generation_time = generation_time.unstack([1, 2]) # use fill_value=... for missing data
    generation_time.columns = new_columns
    generation_time.to_csv("{}_generation_time.csv".format(experiment), na_rep="NaN")
    generation_time.dropna().to_csv("{}_generation_time_clean.csv".format(experiment))

    mean_generation_time = time_accumulated(generation_time)
    mean_generation_time.to_csv("{}_generation_time_accumulated.csv".format(experiment))

    generation_time_renamed = generation_time.rename(mapper=lambda name: name + "_time", axis='columns')
    nodes_count_renamed = nodes_count.rename(mapper=lambda name: name + "_nodes", axis='columns')

    nodes_and_time = pd.concat([generation_time_renamed, nodes_count_renamed], axis=1)
    nodes_and_time.to_csv("{}_generation_all.csv".format(experiment), na_rep="NaN")
    nodes_and_time.dropna().to_csv("{}_generation_all_clean.csv".format(experiment))


if __name__ == '__main__':

    experiments = [("random", 100), ("pldi", 25)]

    speedup_reference = "algorithm0e"

    execution_time_all = []
    gen_results_all = []
    for experiment, n in experiments:

        # execution
        execution_time = read_results_execution(experiment, n, [0, 3]) # [0, 1] for stddev 
        execution_time_all.append(execution_time)
        process_data_execution(execution_time, experiment)

        # generation
        gen_results = read_results_generation(experiment, n)
        gen_results_all.append(gen_results)
        process_data_generation(gen_results, experiment)

    # execution combined
    execution_time_combined = pd.concat(execution_time_all)
    process_data_execution(execution_time_combined, "combined")

    # generation combined
    gen_results_combined = pd.concat(gen_results_all)
    process_data_generation(gen_results_combined, "combined")
 