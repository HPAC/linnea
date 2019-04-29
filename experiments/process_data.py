import os.path
import pandas as pd
import math
from statistics import median

def read_results_generation(experiment, number_of_experiments):

    dirs = []
    new_columns = []
    for dir, col in zip(["c_m", "c_nm", "e_m", "e_nm"], ["constructive_merging", "constructive_no_merging", "exhaustive_merging", "exhaustive_no_merging"]):
        path = os.path.join(experiment, "generation", dir)
        if os.path.exists(path):
            dirs.append(path)
            new_columns.append(col)
        else:
            print("No directory", path)

    file_dfs = []
    for dir in dirs:
        for jobindex in range(1, number_of_experiments+1):
            file_name = "{}/generation{:03}.csv".format(dir, jobindex)
            if os.path.isfile(file_name):
                try:
                    df = pd.read_csv(file_name, index_col=[0, 1, 2])
                except pd.errors.EmptyDataError:
                    print("Empty file", file_name)
                else:
                    file_dfs.append(df)
                
            else:
                print("Missing file", file_name)

    if file_dfs:
        return pd.concat(file_dfs)
    else:
        return None


def ci_pos(n):
    # z = 1.96 # 95%
    z = 2.576 # 99%
    # z = 3.291 # 99.9%
    lower_pos = math.floor((n - z*math.sqrt(n))/2)
    upper_pos = math.ceil(1 + (n + z*math.sqrt(n))/2)
    return lower_pos-1, upper_pos-1


def read_results_execution(experiment_name, number_of_experiments):

    # TODO move this into the loop below
    example_names = ["{}{:03}".format(experiment_name, i) for i in range(1, number_of_experiments+1)]

    language_dfs = []
    for language in ["julia", "matlab", "cpp"]:
        # new
        path = os.path.join(experiment_name, "execution", language)
        if os.path.exists(path):
            file_dfs = []
            for example in example_names:
                file_path = "{0}/{1}_results_{2}_timings.txt".format(path, language, example)
                if os.path.isfile(file_path):
                    # TODO what if rows have different lengths?
                    if language == "matlab":
                        skiprows = 1
                    else:
                        skiprows = 0

                    try:
                        df = pd.read_csv(file_path, sep='\t', skipinitialspace=True, skiprows=skiprows, header=None, index_col=0)
                    except pd.errors.EmptyDataError:
                        print("Empty file", file_path)
                    else:
                        df_processed = pd.DataFrame()
                        df_processed["min_time"] = df.min(axis=1)
                        df_processed["median_time"] = df.median(axis=1)
                        df_processed["mean_time"] = df.mean(axis=1)
                        df_processed["ci_lower"] = df.apply(lambda row: sorted(row)[ci_pos(len(row))[0]], axis=1)
                        df_processed["ci_upper"] = df.apply(lambda row: sorted(row)[ci_pos(len(row))[1]], axis=1)
                        df_processed["ci_lower_rel"] = df_processed["ci_lower"]/df_processed["median_time"]
                        df_processed["ci_upper_rel"] = df_processed["ci_upper"]/df_processed["median_time"]
                        df_processed["min_time_rel"] = df_processed["min_time"]/df_processed["median_time"]
                        df_processed.index = pd.MultiIndex.from_product([[example], df.index.values], names=['example', 'implementation'])
                        file_dfs.append(df_processed)
                else:
                    print("Missing file", file_path)

            df = pd.concat(file_dfs, sort=True)
            language_dfs.append(df)
        else:
            print("No results for", language)

    return pd.concat(language_dfs, sort=True)


def read_intensity(experiment_name, number_of_experiments):
    example_names = ["{}{:03}".format(experiment_name, i) for i in range(1, number_of_experiments+1)]

    file_dfs = []
    for strategy in ["c", "e"]:
        for example_name in example_names:
            file_path = os.path.join(experiment_name, "intensity", strategy, "{}_intensity.csv".format(example_name))
            if os.path.isfile(file_path):
                try:
                    df = pd.read_csv(file_path, index_col=[0, 1])
                except pd.errors.EmptyDataError:
                    print("Empty file", file_path)
                else:
                    file_dfs.append(df)
                
            else:
                print("Missing file", file_path)

    return pd.concat(file_dfs, sort=True)


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

    armadillo = time_data.apply(lambda row: select_best(row["naive_armadillo"], row["recommended_armadillo"]), axis=1)
    matlab = time_data.apply(lambda row: select_best(row["naive_matlab"], row["recommended_matlab"]), axis=1)
    julia = time_data.apply(lambda row: select_best(row["naive_julia"], row["recommended_julia"]), axis=1)
    eigen = time_data.apply(lambda row: select_best(row["naive_eigen"], row["recommended_eigen"]), axis=1)

    res = pd.concat([time_data["algorithm0e"], armadillo, matlab, julia, eigen], axis=1)
    res.rename(columns={0: "armadillo", 1: "matlab", 2: "julia", 3: "eigen"}, inplace=True)

    return res


def to_speedup_data(time_data, reference):
    speedup_data = time_data.apply(lambda row: row/row[reference], axis=1)
    speedup_data["mean"] = speedup_data.apply(lambda row: row.mean(), axis=1)
    return speedup_data


def speedup_CI(time_data, implementation, reference):
    imp_ci_lower = time_data["{}_ci_lower".format(implementation)]
    imp_ci_upper = time_data["{}_ci_upper".format(implementation)]
    ref_ci_lower = time_data["{}_ci_lower".format(reference)]
    ref_ci_upper = time_data["{}_ci_upper".format(reference)]

    if max(imp_ci_lower, ref_ci_lower) <= min(imp_ci_upper, ref_ci_upper):
        # if CIs overlap, return 1
        return 1.0
    else:
        return time_data["{}_median_time".format(implementation)]/time_data["{}_median_time".format(reference)]

def to_speedup_data_CI(time_data):
    df = pd.DataFrame()
    reference = "algorithm0e"
    for implementation in ["naive_julia", "recommended_julia", "naive_matlab", "recommended_matlab", "naive_eigen", "recommended_eigen", "naive_armadillo", "recommended_armadillo", "algorithm0c"]:
        df[implementation] = time_data.apply(speedup_CI, axis=1, args=(implementation, reference))
    df["intensity_e"] = time_data["intensity_e"]
    df["intensity_c"] = time_data["intensity_c"]
    return df

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


def process_data_execution(execution_results, experiment, intensity_cols):

    execution_results = execution_results.unstack()

    execution_time = execution_results["median_time"]

    execution_results = execution_results.reorder_levels([1, 0], axis=1)
    execution_results.columns = ["_".join(col).strip() for col in execution_results.columns.values]
    execution_results = pd.concat([execution_results, intensity_cols], axis=1, sort=True)
    execution_results.to_csv("{}_execution_results.csv".format(experiment), na_rep="NaN")
    execution_results.dropna().to_csv("{}_execution_results_clean.csv".format(experiment), na_rep="NaN")
    to_speedup_data_CI(execution_results.dropna()).to_csv("{}_speedup_CI.csv".format(experiment), na_rep="NaN")

    execution_time_with_intensity = pd.concat([execution_time, intensity_cols], axis=1, sort=True)
    execution_time_with_intensity.to_csv("{}_execution_time.csv".format(experiment), na_rep="NaN")

    execution_time.dropna().to_csv("{}_execution_time_clean.csv".format(experiment))

    performance_profiles_data = to_performance_profiles_data(performance_profiles_data_reduce(execution_time))
    performance_profiles_data.to_csv("{}_performance_profile.csv".format(experiment), na_rep="NaN")
    
    speedup_data = to_speedup_data(execution_time, speedup_reference)
    speedup_data = pd.concat([speedup_data, intensity_cols], axis=1, sort=True)
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

def process_data_intensity(intensity_data, experiment):

    intensity_data.to_csv("{}_intensity_all.csv".format(experiment), na_rep="NaN")

    intensity_only = intensity_data.unstack(level=1).loc[:, ("intensity")]
    intensity_only.to_csv("{}_intensity.csv".format(experiment), na_rep="NaN")

    intensity_only_sorted = intensity_only.sort_values(by=["algorithm0e", "algorithm0c"])
    intensity_only_sorted.to_csv("{}_intensity_sorted.csv".format(experiment), na_rep="NaN")

    intensity_optimal = pd.DataFrame()
    intensity_optimal["intensity_e"] = intensity_only["algorithm0e"]
    intensity_optimal["intensity_c"] = intensity_only["algorithm0c"]

    return intensity_optimal, intensity_only

if __name__ == '__main__':

    experiments = [("random", 100), ("pldi", 25)]

    speedup_reference = "algorithm0e"

    execution_time_all = []
    gen_results_all = []
    intensity_data_all = []
    execution_results_all = []
    for experiment, n in experiments:

        # intensity
        intensity_data = read_intensity(experiment, n)
        intensity_data_all.append(intensity_data)
        intensity_cols, intensity_k_best = process_data_intensity(intensity_data, experiment)

        # execution
        execution_results = read_results_execution(experiment, n)
        execution_results_all.append(execution_results)
        process_data_execution(execution_results, experiment, intensity_cols)

        # generation
        gen_results = read_results_generation(experiment, n)
        gen_results_all.append(gen_results)
        process_data_generation(gen_results, experiment)

    # intensity combined
    intensity_data_combined = pd.concat(intensity_data_all)
    intensity_cols, intensity_k_best = process_data_intensity(intensity_data_combined, "combined")

    # execution combined
    execution_results_combined = pd.concat(execution_results_all, sort=True)
    process_data_execution(execution_results_combined, "combined", intensity_cols)

    # generation combined
    gen_results_combined = pd.concat(gen_results_all)
    process_data_generation(gen_results_combined, "combined")


 