import os.path
import argparse
import pandas as pd

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

    return final


def to_speedup_data(time_data, reference):
    speedup_data = time_data.apply(lambda row: row/row[reference], axis=1)
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


def read_results(experiment_name, number_of_experiments, usecols=range(0, 2)):

    """ TODO idea: Instead of usecols, use Enum an argument, which is maped to
    columns. This can then also be used in df.rename()

    What if we want more than one column?
    """

    # TODO move this into the loop below
    example_names = ["{}{}".format(experiment_name, i) for i in range(1, number_of_experiments+1)]

    language_dfs = []
    for language in ["cpp", "julia", "matlab"]:
        if os.path.exists(language):
            file_dfs = []
            for example in example_names:
                file_path = "{0}/{0}_results_{1}.txt".format(language, example)
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


def check_std_dev(experiment_name, number_of_experiments, threshold=0.1):
    execution_time = read_results(experiment_name, number_of_experiments)
    std_dev = read_results(experiment_name, number_of_experiments, usecols=[0, 2])
    rel_std_dev = std_dev/execution_time.values
    rel_std_dev = rel_std_dev.stack()
    rel_std_dev = rel_std_dev.loc[rel_std_dev > threshold]
    rel_std_dev.sort_values(ascending=False, inplace=True)
    return rel_std_dev


if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog="process_data")
    parser.add_argument("experiment", choices=["lamp_example", "random"])
    args = parser.parse_args()

    if args.experiment == "lamp_example":
        number_of_experiments = 31
    elif args.experiment == "random":
        number_of_experiments = 100

    # execution_time = read_results()
    execution_time = read_results(args.experiment, number_of_experiments, [0, 3]) # use min
    
    execution_time.to_csv("execution_time.csv", na_rep="NaN")

    performance_profiles_data = to_performance_profiles_data(execution_time)
    # performance_profiles_data.rename(columns={0:"y", "generated0":"linnea"}, inplace=True)
    performance_profiles_data.to_csv("performance_profile.csv", na_rep="NaN")

    speedup_data = to_speedup_data(execution_time, "algorithm0c")
    speedup_data.to_csv("speedup.csv", na_rep="NaN")

    mean_speedup = to_mean_speedup(speedup_data)
    mean_speedup.to_csv("mean_speedup.csv", na_rep="NaN")

    # print(mean_speedup)

    speedup_over_linnea = compare_to_fastest(execution_time, "algorithm0c")
    speedup_over_linnea.to_csv("speedup_over_linnea.csv", na_rep="NaN")
    # print(check_std_dev())
    

