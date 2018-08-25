import pkg_resources
import json
import os.path
import itertools

# this causes a problem because it tries to read the config.json in the current directory.
# from linnea.config import Strategy

def load_config():

    config_file = "jobscripts/config.json"
    if os.path.exists(config_file):
        with open(config_file) as jsonfile:
            data = json.load(jsonfile)

            data_time = dict()
            for key, value in data["time"].items():
                if key == "exclusive":
                    if value:
                        data_time["exclusive"] = "#BSUB -x                     # exclusive access"
                    else:
                        data_time["exclusive"] = ""
                else:
                    data_time[key] = value

    return data["generate"], data_time

def time_generation_script(replacement):

    template_path = "jobscripts/templates/time_generation.sh"
    template_str = pkg_resources.resource_string(__name__, template_path).decode("UTF-8")

    for strategy, merging in itertools.product(["c", "e"], [True, False]):
        replacement_copy = replacement.copy()

        file_name_parts = ["time_generation"]

        file_name_parts.append(strategy)
        replacement_copy["strategy"] = strategy
        # TODO there is some redundancy here
        if merging:
            replacement_copy["strategy_name"] = strategy + "_m"
        else:
            replacement_copy["strategy_name"] = strategy + "_nm"
        # if strategy is Strategy.constructive:
        #     file_name_parts.append("c")
        #     replacement_copy["strategy"] = "c"
        # elif strategy is Strategy.exhaustive:
        #     file_name_parts.append("e")
        #     replacement_copy["strategy"] = "e"

        if merging:
            file_name_parts.append("m")
            replacement_copy["merging"] = "true"
        else:
            file_name_parts.append("nm")
            replacement_copy["merging"] = "false"

        file_name = "jobscripts/{}/{}.sh".format(replacement_copy["name"], "_".join(file_name_parts))
        with open(file_name, "wt", encoding='utf-8') as output_file:
            print("Writing", file_name)
            output_file.write(template_str.format(**replacement_copy))


def time_execution_scripts(replacement):

    for language in ["julia", "matlab", "cpp"]:
        template_path = "jobscripts/templates/time_{}.sh".format(language)
        template_str = pkg_resources.resource_string(__name__, template_path).decode("UTF-8")

        file_name = "jobscripts/{}/time_{}.sh".format(replacement["name"], language)
        with open(file_name, "wt", encoding='utf-8') as output_file:
            print("Writing", file_name)
            output_file.write(template_str.format(**replacement))


def generate_code_scripts(replacement):

    template_path = "jobscripts/templates/generate_code.sh"
    template_str = pkg_resources.resource_string(__name__, template_path).decode("UTF-8")

    for strategy in ["c", "e"]:
        replacement_copy = replacement.copy()

        replacement_copy["strategy"] = strategy

        file_name = "jobscripts/{}/generate_code_{}.sh".format(replacement_copy["name"], strategy)
        with open(file_name, "wt", encoding='utf-8') as output_file:
            print("Writing", file_name)
            output_file.write(template_str.format(**replacement_copy))


def generate_scripts(experiment, number_of_experiments):

    replacement_generate, replacement_time = load_config()

    replacement_generate["jobs"] = number_of_experiments
    replacement_time["jobs"] = number_of_experiments
    replacement_generate["name"] = experiment
    replacement_time["name"] = experiment

    time_generation_script(replacement_time)
    time_execution_scripts(replacement_time)
    generate_code_scripts(replacement_generate)


if __name__ == '__main__':
    generate_scripts()