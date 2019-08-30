import os.path
import pandas as pd
import math
from statistics import median

def read_results_generation(experiment, number_of_experiments):

    dirs = []
    new_columns = []
    for merging in "merging", "no_merging":
        path = os.path.join(experiment, "generation", merging)
        if os.path.exists(path):
            dirs.append(path)
            new_columns.append(merging)
        else:
            print("No directory", path)

    file_dfs = []
    for dir in dirs:
        for jobindex in range(1, number_of_experiments+1):
            file_name = "{}/generation{:03}.csv".format(dir, jobindex)
            if os.path.isfile(file_name):
                try:
                    df = pd.read_csv(file_name, index_col=[0, 1])
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


def read_execution_file(file, skiprows):
    """Function to read data with variable number of columns.

    """
    from io import StringIO
    from os import stat
    if os.stat(file).st_size != 0:
        with open(file) as input_file:
            file_dfs = []
            for i, line in enumerate(input_file):
                if i < skiprows:
                    continue
                df = pd.read_csv(StringIO(line), sep='\t', skipinitialspace=True, header=None, index_col=[0, 1])
                file_dfs.append(df)
            df = pd.concat(file_dfs)
            # print(df)
        return df 
    else:
        # print("Empty file", file)
        return None


def read_results_execution(experiment_name, number_of_experiments):

    # TODO move this into the loop below
    example_names = ["{}{:03}".format(experiment_name, i) for i in range(1, number_of_experiments+1)]

    language_dfs = []
    for language in ["julia", "matlab", "cpp"]:
        for threads in [1, 24]:
            path = os.path.join(experiment_name, "execution", language, "t{}".format(threads))
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

                        # df = read_execution_file(file_path, skiprows)
                        # if df is None:
                        #     print("Empty file", file_path)
                        try:
                            df = pd.read_csv(file_path, sep='\t', skipinitialspace=True, skiprows=skiprows, header=None, index_col=[0, 1])
                        except pd.errors.EmptyDataError:
                            print("Empty file", file_path)
                        else:
                            df_processed = pd.DataFrame()
                            df_processed["min_time"] = df.min(axis=1)
                            df_processed["median_time"] = df.median(axis=1)
                            df_processed["mean_time"] = df.mean(axis=1)
                            df_processed["ci_lower"] = df.apply(lambda row: sorted(row.dropna())[ci_pos(len(row.dropna()))[0]], axis=1)
                            df_processed["ci_upper"] = df.apply(lambda row: sorted(row.dropna())[ci_pos(len(row.dropna()))[1]], axis=1)
                            df_processed["ci_lower_rel"] = df_processed["ci_lower"]/df_processed["median_time"]
                            df_processed["ci_upper_rel"] = df_processed["ci_upper"]/df_processed["median_time"]
                            df_processed["min_time_rel"] = df_processed["min_time"]/df_processed["median_time"]
                            df_processed["iters"] = df.apply(lambda row: len(row.dropna()), axis=1)
                            tuples = [(example, implementation, threads) for implementation, threads in df.index.values]
                            df_processed.index = pd.MultiIndex.from_tuples(tuples, names=['example', 'implementation', 'threads'])
                            file_dfs.append(df_processed)
                    else:
                        print("Missing file", file_path)

                if file_dfs:
                    df = pd.concat(file_dfs, sort=True)
                    language_dfs.append(df)
            else:
                print("No directory", path)

    if language_dfs:
        return pd.concat(language_dfs, sort=True)
    else:
        return pd.DataFrame()


def read_intensity(experiment_name, number_of_experiments, dir):
    example_names = ["{}{:03}".format(experiment_name, i) for i in range(1, number_of_experiments+1)]

    file_dfs = []
    for example_name in example_names:
        file_path = os.path.join(experiment_name, dir, "{}_intensity.csv".format(example_name))
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


def read_trace(experiment_name, number_of_experiments):
    example_names = ["{}{:03}".format(experiment_name, i) for i in range(1, number_of_experiments+1)]
    file_dfs = []
    for merging, merging_bool in [("merging", True), ("no_merging", False)]:
        path = file_path = os.path.join(experiment_name, "generation", merging)
        if os.path.exists(path) and os.listdir(path):
            for example in example_names:
                file_path = os.path.join(path, "{}_trace.csv".format(example))
                if os.path.isfile(file_path):
                    try:
                        df = pd.read_csv(file_path, index_col=0)
                    except pd.errors.EmptyDataError:
                        print("Empty file", file_path)
                    else:
                        processed_data = dict()
                        # processed_data = pd.DataFrame() # [], index=[], columns=[example]
                        min_cost = df["cost"].min()
                        max_cost = df["cost"].max()
                        processed_data["first_solution_time"] = df["time"].min()
                        processed_data["best_solution_time"] = df["time"].max()
                        processed_data["first_solution_cost"] = max_cost
                        processed_data["best_solution_cost"] = min_cost
                        for p in [1, 5, 10, 15, 25]:
                            processed_data["time_to_{}pc".format(p)] = df.loc[df["cost"] < (min_cost * (1+p/100))]["time"].min()

                        mindex = pd.MultiIndex.from_tuples([(example, merging_bool)], names=['example', 'merging'])

                        file_dfs.append(pd.DataFrame(processed_data, index=mindex))
                else:
                    print("Missing file", file_path)
        else:
            print(path, "is empty or does not exist.")

    if file_dfs:
        return pd.concat(file_dfs, sort=True)
    else:
        return pd.DataFrame()


def read_write_time(experiment_name, number_of_experiments):
    example_names = ["{}{:03}".format(experiment_name, i) for i in range(1, number_of_experiments+1)]
    file_dfs = []
    for merging, merging_bool in [("merging", True), ("no_merging", False)]:
        path = file_path = os.path.join(experiment_name, "generation", merging)
        if os.path.exists(path) and os.listdir(path):
            for example in example_names:
                file_path = os.path.join(path, "{}_code_gen_time.csv".format(example))
                if os.path.isfile(file_path):
                    try:
                        df = pd.read_csv(file_path, index_col=[0, 1])
                    except pd.errors.EmptyDataError:
                        print("Empty file", file_path)
                    else:
                        file_dfs.append(df)
                else:
                    print("Missing file", file_path)
        else:
            print(path, "is empty or does not exist.")

    if file_dfs:
        return pd.concat(file_dfs, sort=True)
    else:
        return pd.DataFrame()


def read_results_k_best(experiment_name, number_of_experiments):

    example_names = ["{}{:03}".format(experiment_name, i) for i in range(1, number_of_experiments+1)]

    path = os.path.join(experiment_name, "execution", "k_best")
    file_dfs = []
    file_dfs_raw = dict()
    for threads in [1, 24]:
        file_dfs_raw[threads] = []

    for example in example_names:
        for threads in [1, 24]:

            file_path = "{0}/t{1}/julia_results_{2}_timings.txt".format(path, threads, example)
            if os.path.isfile(file_path):
                try:
                    df = pd.read_csv(file_path, sep='\t', skipinitialspace=True, header=None, index_col=[0, 1])
                except pd.errors.EmptyDataError:
                    print("Empty file", file_path)
                else:
                    # TODO put this into a function
                    df_processed = pd.DataFrame()
                    df_processed["min_time"] = df.min(axis=1)
                    df_processed["median_time"] = df.median(axis=1)
                    df_processed["mean_time"] = df.mean(axis=1)
                    df_processed["ci_lower"] = df.apply(lambda row: sorted(row)[ci_pos(len(row))[0]], axis=1)
                    df_processed["ci_upper"] = df.apply(lambda row: sorted(row)[ci_pos(len(row))[1]], axis=1)
                    df_processed["ci_lower_rel"] = df_processed["ci_lower"]/df_processed["median_time"]
                    df_processed["ci_upper_rel"] = df_processed["ci_upper"]/df_processed["median_time"]
                    df_processed["min_time_rel"] = df_processed["min_time"]/df_processed["median_time"]
                    tuples = [(example, implementation, threads) for implementation, threads in df.index.values]
                    df_processed.index = pd.MultiIndex.from_tuples(tuples, names=['example', 'implementation', 'threads'])
                    file_dfs.append(df_processed)

                    file_dfs_raw[threads].append(pd.DataFrame(df.drop(["naive_julia", "recommended_julia"], level=0, errors='ignore').stack().reset_index(drop=True), columns=[example]))
            else:
                print("Missing file", file_path)

    if file_dfs:
        for threads, data in file_dfs_raw.items():
            df_raw = pd.concat(data, axis=1)
            df_raw.to_csv("{}_t{}_k_best_raw.csv".format(experiment, threads), na_rep="NaN")

        df = pd.concat(file_dfs, sort=True)
        return df
    else:
        print("No results for k best experiment.")
        return pd.DataFrame()


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

    res = pd.concat([time_data["algorithm0"], armadillo, matlab, julia, eigen], axis=1)
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
    reference = "algorithm0"
    for implementation in ["naive_julia", "recommended_julia", "naive_matlab", "recommended_matlab", "naive_eigen", "recommended_eigen", "naive_armadillo", "recommended_armadillo"]:
        df[implementation] = time_data.apply(speedup_CI, axis=1, args=(implementation, reference))
    df["intensity"] = time_data["intensity"]
    return df


def to_mean_speedup(speedup_data, drop_reference=None):
    speedup_mean = pd.DataFrame(speedup_data.mean())
    if drop_reference:
        speedup_mean.drop(drop_reference, inplace=True)
    speedup_mean = speedup_mean.transpose()
    return speedup_mean


def compare_to_fastest(speedup_data, intensity_data):
    """Identifies cases where other systems are faster than reference.

    Returns a DataFrame with example problem, algorithm and slowdown, sorted by
    the magintude of the slowdown.
    """
    df = speedup_data.drop(["algorithm0", "mean", "intensity"], axis=1, errors="ignore")
    df = 1. / df 
    df = df.stack()
    df = pd.DataFrame(df, columns=["speedup"]).join(intensity_data, how='inner')
    df = df.loc[df["speedup"] > 1]
    df.sort_values(by=['speedup'], ascending=False, inplace=True)
    return df


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
        experiments = ["random", "application"]
    else:    
        experiments = [experiment]

    # TODO this code also appears in "read_results_generation"    
    new_columns = []
    for exp in experiments:
        for merging in ["merging", "no_merging"]:
            if merging in new_columns:
                continue
            path = os.path.join(exp, "generation", merging)
            if os.path.exists(path):
                new_columns.append(merging)

    return new_columns


def process_data_execution(data, experiment_name, intensity_cols):

    # for threads in [24]:
    for threads in [1, 24]:
        experiment = "{}_t{}".format(experiment_name, threads)

        execution_results = data.xs(threads, level=2)
        execution_results = execution_results.unstack()

        execution_time = execution_results[base_time]

        execution_results = execution_results.reorder_levels([1, 0], axis=1)
        execution_results.columns = ["_".join(col).strip() for col in execution_results.columns.values]
        execution_results = pd.concat([execution_results, intensity_cols], axis=1, sort=True)
        execution_results.to_csv("{}_execution_results.csv".format(experiment), na_rep="NaN")
        execution_results.dropna().to_csv("{}_execution_results_clean.csv".format(experiment), na_rep="NaN")
        speedup_data_CI = to_speedup_data_CI(execution_results.dropna())
        speedup_data_CI.to_csv("{}_speedup_CI.csv".format(experiment), na_rep="NaN")

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

        speedup_over_linnea = compare_to_fastest(speedup_data, intensity_cols)
        speedup_over_linnea.to_csv("{}_speedup_over_linnea.csv".format(experiment), na_rep="NaN")

        speedup_over_linnea_CI = compare_to_fastest(speedup_data_CI, intensity_cols)
        speedup_over_linnea_CI.to_csv("{}_speedup_over_linnea_CI.csv".format(experiment), na_rep="NaN")


def process_data_generation(gen_results, experiment):

    gen_results.to_csv("{}_generation.csv".format(experiment), na_rep="NaN")

    new_columns = get_column_names(experiment)

    nodes_count = gen_results.drop(columns=["mean", "std", "min", "max", "solution"])
    nodes_count = nodes_count.unstack(1) # use fill_value=... for missing data
    nodes_count.columns = new_columns
    nodes_count.to_csv("{}_generation_nodes_count.csv".format(experiment), na_rep="NaN")
    nodes_count.dropna().to_csv("{}_generation_nodes_count_clean.csv".format(experiment))

    generation_time = gen_results.drop(columns=["std", "min", "max", "nodes", "solution"])
    generation_time = generation_time.unstack(1) # use fill_value=... for missing data
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


def process_data_intensity(intensity_data, experiment, output_name):

    intensity_data.to_csv("{}_{}_all.csv".format(experiment, output_name), na_rep="NaN")

    intensity_only = intensity_data.unstack(level=1).loc[:, ("intensity")]
    intensity_only.to_csv("{}_{}.csv".format(experiment, output_name), na_rep="NaN")

    # what is the point of this?
    # intensity_only_sorted = intensity_only.sort_values(by=["algorithm0"])
    # intensity_only_sorted.to_csv("{}_{}_sorted.csv".format(experiment, output_name), na_rep="NaN")

    intensity_optimal = pd.DataFrame()
    intensity_optimal["intensity"] = intensity_only["algorithm0"]

    return intensity_optimal


def process_data_trace(trace_data, experiment):
    if not trace_data.empty:
        trace_data.to_csv("{}_generation.csv".format(experiment), na_rep="NaN")

        for merging, merging_bool in [("merging", True), ("no_merging", False)]:
            try:
                data = trace_data.xs(merging_bool, level=1)
            except KeyError:
                pass
            else:
                data.to_csv("{}_generation_{}.csv".format(experiment, merging), na_rep="NaN")


def process_write_time(write_time_data, experiment):
    if not write_time_data.empty:
        write_time_data.to_csv("{}_write_time.csv".format(experiment), na_rep="NaN")
        write_time_data.unstack().mean().to_csv("{}_write_time_avg.csv".format(experiment), na_rep="NaN", header=True)


def process_data_k_best(data, intensity_data, experiment_name):

    for threads in [1, 24]:
        experiment = "{}_t{}".format(experiment_name, threads)

        k_best_data = data.xs(threads, level=2)

        plain_data = pd.concat([k_best_data, intensity_data], axis=1, sort=True)
        plain_data.drop(["naive_julia", "recommended_julia"], level=1, inplace=True)
        plain_data.to_csv("{}_k_best.csv".format(experiment), na_rep="NaN")

        for col in ["min_time", "min_time_rel", "median_time", "ci_upper_rel", "ci_lower_rel", "cost"]:
            df = plain_data[col].unstack(0)
            df.index = df.index.map(lambda x: int(x[9:]))
            df.sort_index(inplace=True)
            df.to_csv("{}_k_best_{}.csv".format(experiment, col), na_rep="NaN")

        # stats algorithm0
        algorithm0_stats = pd.DataFrame(plain_data.xs("algorithm0", level=1))
        algorithm0_stats.dropna(inplace=True)

        df = plain_data["min_time"].unstack(0)
        algorithm0_stats["k"] = df.count()

        df = plain_data["cost"].unstack(0)
        algorithm0_stats["cost_cutoff"] = df.max()/df.min()
        
        algorithm0_stats.to_csv("{}_k_best_stats_algorithm0.csv".format(experiment), na_rep="NaN")

        # stats for solution with smallest min time
        min_time_stats = plain_data.loc[plain_data.groupby(level=0).min_time.idxmin().dropna()]
        min_time_stats.reset_index(level=1, inplace=True)

        min_time_stats["idx"] = min_time_stats.apply(lambda row: int(row["implementation"][9:]), axis=1)

        del min_time_stats["ci_lower"]
        del min_time_stats["ci_lower_rel"]
        del min_time_stats["ci_upper"]
        del min_time_stats["ci_upper_rel"]
        del min_time_stats["min_time_rel"]

        min_time_stats["slowdown"] = min_time_stats["min_time"]/algorithm0_stats["min_time"]
        min_time_stats["speedup"] = algorithm0_stats["min_time"]/min_time_stats["min_time"]
        min_time_stats["cost_rel"] = min_time_stats["cost"]/algorithm0_stats["cost"]
        min_time_stats["intensity_ref"] = algorithm0_stats["intensity"]

        min_time_stats.to_csv("{}_k_best_stats_min_time.csv".format(experiment), na_rep="NaN")

        # stats for solution with smallest median time
        median_time_stats = plain_data.loc[plain_data.groupby(level=0).median_time.idxmin().dropna()]
        median_time_stats.reset_index(level=1, inplace=True)
        median_time_stats["idx"] = median_time_stats.apply(lambda row: int(row["implementation"][9:]), axis=1)

        median_time_stats["slowdown"] = median_time_stats["median_time"]/algorithm0_stats["median_time"]
        median_time_stats["speedup"] = algorithm0_stats["median_time"]/median_time_stats["median_time"]
        tmp = pd.DataFrame()
        tmp["optimal_ci_lower"] = median_time_stats["ci_lower"]
        tmp["optimal_ci_upper"] = median_time_stats["ci_upper"]
        tmp["optimal_median_time"] = median_time_stats["median_time"]
        tmp["algorithm0_ci_lower"] = algorithm0_stats["ci_lower"]
        tmp["algorithm0_ci_upper"] = algorithm0_stats["ci_upper"]
        tmp["algorithm0_median_time"] = algorithm0_stats["median_time"]
        median_time_stats["slowdown_CI"] = tmp.apply(speedup_CI, axis=1, args=("optimal", "algorithm0"))
        median_time_stats["speedup_CI"] = tmp.apply(speedup_CI, axis=1, args=("algorithm0", "optimal"))
        median_time_stats["cost_rel"] = median_time_stats["cost"]/algorithm0_stats["cost"]
        median_time_stats["intensity_ref"] = algorithm0_stats["intensity"]

        median_time_stats.to_csv("{}_k_best_stats_median_time.csv".format(experiment), na_rep="NaN")

        # TODO search for cases where speedup < 1  and cost_rel > 1
        # median_time_stats.loc[(median_time_stats["speedup"] < 1.0) & (median_time_stats["cost_rel"] > 1.0)]


def process_data_optimal_solution(execution_data, k_best_data, intensity_data, experiment_name):

    # for threads in [24]:
    for threads in [1, 24]:
        experiment = "{}_t{}".format(experiment_name, threads)

        execution_results = execution_data.xs(threads, level=2)
        execution_results = execution_results.unstack()

        execution_time = execution_results[base_time].copy()

        execution_results = execution_results.reorder_levels([1, 0], axis=1)
        execution_results.columns = ["_".join(col).strip() for col in execution_results.columns.values]
        execution_results = pd.concat([execution_results, intensity_cols], axis=1, sort=True)

        k_best = k_best_data.xs(threads, level=2).copy()
        k_best.drop(["naive_julia", "recommended_julia"], level=1, inplace=True)

        execution_time["algorithm0"] = k_best.xs("algorithm0", level=1)[base_time]

        speedup_data = to_speedup_data(execution_time, speedup_reference)
        speedup_over_linnea = compare_to_fastest(speedup_data, intensity_cols)
        speedup_over_linnea.to_csv("{}_k_best_speedup_over_linnea.csv".format(experiment), na_rep="NaN")

        k_best_time_stats = k_best.loc[k_best.groupby(level=0)[base_time].idxmin().dropna()]
        k_best_time_stats.reset_index(level=1, inplace=True)

        execution_time["algorithm0"] = k_best_time_stats[base_time]

        speedup_data = to_speedup_data(execution_time, speedup_reference)
        speedup_data = pd.concat([speedup_data, intensity_cols], axis=1, sort=True)
        speedup_data.to_csv("{}_speedup_OS.csv".format(experiment), na_rep="NaN")

        execution_results["algorithm0_min_time"] = k_best_time_stats["min_time"]
        execution_results["algorithm0_mean_time"] = k_best_time_stats["mean_time"]
        execution_results["algorithm0_median_time"] = k_best_time_stats["median_time"]
        execution_results["algorithm0_ci_lower"] = k_best_time_stats["ci_lower"]
        execution_results["algorithm0_ci_upper"] = k_best_time_stats["ci_upper"]

        speedup_data_CI = to_speedup_data_CI(execution_results.dropna())
        speedup_data_CI.to_csv("{}_speedup_OS_CI.csv".format(experiment), na_rep="NaN")

        # execution_time = pd.concat([execution_time, intensity_cols], axis=1, sort=True)
        
        performance_profiles_data = to_performance_profiles_data(performance_profiles_data_reduce(execution_time))
        performance_profiles_data.to_csv("{}_performance_profile_OS.csv".format(experiment), na_rep="NaN")
        
        speedup_over_linnea = compare_to_fastest(speedup_data, intensity_cols)
        speedup_over_linnea.to_csv("{}_k_best_speedup_over_linnea_OS.csv".format(experiment), na_rep="NaN")

        speedup_over_linnea_CI = compare_to_fastest(speedup_data_CI, intensity_cols)
        speedup_over_linnea_CI.to_csv("{}_k_best_speedup_over_linnea_OS_CI.csv".format(experiment), na_rep="NaN")


if __name__ == '__main__':

    experiments = [("random", 100), ("application", 25)]

    base_time = "min_time" # what is used for calculating speedups
    speedup_reference = "algorithm0"

    execution_time_all = []
    gen_results_all = []
    intensity_data_all = []
    intensity_k_best_data_all = []
    trace_data_all = []
    write_data_all = []
    k_best_data_all = []
    execution_results_all = []
    for experiment, n in experiments:

        # trace
        trace_data = read_trace(experiment, n)
        trace_data_all.append(trace_data)
        process_data_trace(trace_data, experiment)

        # write time
        write_data = read_write_time(experiment, n)
        write_data_all.append(write_data)
        process_write_time(write_data, experiment)

        # intensity k best
        intensity_data = read_intensity(experiment, n, "intensity_k_best")
        intensity_k_best_data_all.append(intensity_data)
        intensity_cols = process_data_intensity(intensity_data, experiment, "intensity_k_best")

        # k best
        k_best_data = read_results_k_best(experiment, n)
        k_best_data_all.append(k_best_data)
        if not k_best_data.empty:
            process_data_k_best(k_best_data, intensity_data, experiment)

        # intensity
        intensity_data = read_intensity(experiment, n, "intensity")
        intensity_data_all.append(intensity_data)
        intensity_cols = process_data_intensity(intensity_data, experiment, "intensity")

        # execution
        execution_results = read_results_execution(experiment, n)
        execution_results_all.append(execution_results)
        if not execution_results.empty:
            process_data_execution(execution_results, experiment, intensity_cols)

        # # generation
        # gen_results = read_results_generation(experiment, n)
        # gen_results_all.append(gen_results)
        # process_data_generation(gen_results, experiment)

        # speedup plot and performance profile with optimal solution as reference
        process_data_optimal_solution(execution_results, k_best_data, intensity_cols, experiment)


    # trace combined
    trace_data_combined = pd.concat(trace_data_all, sort=True)
    process_data_trace(trace_data_combined, "combined")

    # write time
    write_data_combined = pd.concat(write_data_all, sort=True)
    process_write_time(write_data_combined, "combined")

    # intensity k best combined
    intensity_data_combined = pd.concat(intensity_k_best_data_all)
    intensity_cols = process_data_intensity(intensity_data_combined, "combined", "intensity_k_best")

    # k best
    k_best_data_combined = pd.concat(k_best_data_all)
    if not k_best_data_combined.empty:
        process_data_k_best(k_best_data_combined, intensity_data_combined, "combined")

    # intensity combined
    intensity_data_combined = pd.concat(intensity_data_all)
    intensity_cols = process_data_intensity(intensity_data_combined, "combined", "intensity")

    # execution combined
    execution_results_combined = pd.concat(execution_results_all, sort=True)
    if not execution_results.empty:
        process_data_execution(execution_results_combined, "combined", intensity_cols)

    # # generation combined
    # gen_results_combined = pd.concat(gen_results_all)
    # process_data_generation(gen_results_combined, "combined")

    # speedup plot and performance profile with optimal solution as reference
    process_data_optimal_solution(execution_results_combined, k_best_data_combined, intensity_cols, "combined")


 