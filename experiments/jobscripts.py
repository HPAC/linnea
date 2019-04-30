import pkg_resources
import os.path
import os
import itertools
from linnea import config

def time_generation_script(replacement):

    template_path = "jobscripts/templates/time_generation.sh"
    template_str = pkg_resources.resource_string(__name__, template_path).decode("UTF-8")

    for strategy, merging in itertools.product(["c", "e"], [True, False]):
        replacement_copy = replacement.copy()

        file_name_parts = ["time_generation"]

        file_name_parts.append(strategy)
        replacement_copy["strategy"] = strategy
        if merging:
            replacement_copy["strategy_name"] = strategy + "_m"
            file_name_parts.append("m")
            replacement_copy["merging"] = "true"
        else:
            replacement_copy["strategy_name"] = strategy + "_nm"
            file_name_parts.append("nm")
            replacement_copy["merging"] = "false"

        file_name = "{}/{}/{}.sh".format(replacement_copy['linnea_jobscripts_path'], replacement_copy["name"],
                                         "_".join(file_name_parts))
        with open(file_name, "wt", encoding='utf-8') as output_file:
            print("Writing", file_name)
            output_file.write(template_str.format(**replacement_copy))


def time_execution_scripts(replacement):

    for language in ["julia", "matlab", "cpp"]:
        template_path = "jobscripts/templates/time_{}.sh".format(language)
        template_str = pkg_resources.resource_string(__name__, template_path).decode("UTF-8")

        file_name = "{}/{}/time_{}.sh".format(replacement['linnea_jobscripts_path'], replacement["name"], language)
        with open(file_name, "wt", encoding='utf-8') as output_file:
            print("Writing", file_name)
            output_file.write(template_str.format(runner_name="runner", output_subdir=language, **replacement))

def time_k_best_script(replacement):

    template_path = "jobscripts/templates/time_julia.sh"
    template_str = pkg_resources.resource_string(__name__, template_path).decode("UTF-8")

    file_name = "{}/{}/time_k_best.sh".format(replacement['linnea_jobscripts_path'], replacement["name"])
    with open(file_name, "wt", encoding='utf-8') as output_file:
        print("Writing", file_name)
        output_file.write(template_str.format(runner_name="runner_k_best", output_subdir="k_best", **replacement))

def generate_code_scripts(replacement):

    template_path = "jobscripts/templates/generate_code.sh"
    template_str = pkg_resources.resource_string(__name__, template_path).decode("UTF-8")

    for arg in ["c", "e", "f"]:
        replacement_copy = replacement.copy()

        replacement_copy["args"] = arg

        if arg == "f":
            replacement_copy["compile"] = "true"
        else:
            replacement_copy["compile"] = "false"

        file_name = "{}/{}/generate_code_{}.sh".format(replacement_copy['linnea_jobscripts_path'],
                                                       replacement_copy["name"], arg)
        with open(file_name, "wt", encoding='utf-8') as output_file:
            print("Writing", file_name)
            output_file.write(template_str.format(**replacement_copy))

scheduler_vars = {
    "LSF":
        {
        "directive":            "BSUB",
        "flag_jobname":         "-J",
        "flag_output":          "-oo",
        "flag_time":            "-W",
        "flag_memory":          "-M ",
        "flag_group":           "-P",
        "flag_model":           "-R",
        "flag_exclusive":       "-x",
        "var_array_idx":        "LSB_JOBINDEX",
        "string_array_idx":     "%I"
        },
    "SLURM":
        {
        "directive":            "SBATCH",
        "flag_jobname":         "-J",
        "flag_output":          "-o",
        "flag_time":            "-t",
        "flag_memory":          "--mem=",
        "flag_group":           "-A",
        "flag_model":           "-p",
        "flag_exclusive":       "--exclusive",
        "var_array_idx":        "SLURM_ARRAY_TASK_ID",
        "string_array_idx":     "%a"
        }
    }

def generate_scripts(experiment, number_of_experiments):

    experiment_configuration = config.experiment_configuration

    scheduler = experiment_configuration["scheduler"]

    if experiment_configuration["time"]["exclusive"]:
        experiment_configuration["time"]["spec_exclusive"] = "#{directive} {flag_exclusive}".format(**scheduler_vars[scheduler])
    else:
        experiment_configuration["time"]["spec_exclusive"] = ""

    for mode in ["time", "generate"]:
        if scheduler == "LSF":
            experiment_configuration[mode]["lsf_arrayjob"] = "[1-{}]".format(number_of_experiments)
            experiment_configuration[mode]["slurm_arrayjob"] = ""
            experiment_configuration[mode]["spec_model"] = "#BSUB -R {model}".format(**experiment_configuration[mode])
        elif scheduler == "SLURM":
            experiment_configuration[mode]["lsf_arrayjob"] = ""
            experiment_configuration[mode]["slurm_arrayjob"] = "#SBATCH --array=1-{}".format(number_of_experiments)
            experiment_configuration[mode]["spec_model"] = ""
        else:
            pass # TODO throw exception?

        experiment_configuration[mode]['name'] = experiment


    dirname = "{}/{}/".format(experiment_configuration['path']['linnea_jobscripts_path'], experiment)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    time_configuration = {**experiment_configuration['time'], **experiment_configuration['path'], **experiment_configuration['version'], **scheduler_vars[scheduler]}
    generate_configuration = {**experiment_configuration['generate'], **experiment_configuration['path'], **experiment_configuration['version'], **scheduler_vars[scheduler]}

    time_generation_script(time_configuration)
    time_execution_scripts(time_configuration)
    time_k_best_script(time_configuration)
    generate_code_scripts(generate_configuration)

if __name__ == '__main__':
    generate_scripts()
