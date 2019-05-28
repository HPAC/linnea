import pkg_resources
import os.path
import os
import itertools
from linnea import config

def time_generation_script(replacement):

    template_path = "jobscripts/templates/time_generation.sh"
    template_str = pkg_resources.resource_string(__name__, template_path).decode("UTF-8")

    for merging in [True, False]:
        replacement_copy = replacement.copy()

        file_name_parts = ["time_generation"]

        if merging:
            replacement_copy["dir_name"] = "merging"
            replacement_copy["merging_label"] = "m"
            file_name_parts.append("m")
            replacement_copy["merging"] = "true"
        else:
            replacement_copy["dir_name"] = "no_merging"
            replacement_copy["merging_label"] = "nm"
            file_name_parts.append("nm")
            replacement_copy["merging"] = "false"

        file_name = "{}/{}/time_generation_{}.sh".format(
                                replacement_copy['linnea_jobscripts_path'],
                                replacement_copy["name"],
                                replacement_copy["merging_label"]
                            )
        with open(file_name, "wt", encoding='utf-8') as output_file:
            print("Writing", file_name)
            output_file.write(template_str.format(**replacement_copy))


def time_execution_scripts(replacement):

    for language in ["julia", "matlab", "cpp"]:
        template_path = "jobscripts/templates/time_{}.sh".format(language)
        template_str = pkg_resources.resource_string(__name__, template_path).decode("UTF-8")

        file_name = "{}/{}/time_{}_t{}.sh".format(replacement['linnea_jobscripts_path'], replacement["name"], language, replacement["threads"])
        with open(file_name, "wt", encoding='utf-8') as output_file:
            print("Writing", file_name)
            output_file.write(template_str.format(runner_name="runner_t{}".format(replacement["threads"]), output_subdir=language, **replacement))

def time_k_best_script(replacement):

    template_path = "jobscripts/templates/time_julia.sh"
    template_str = pkg_resources.resource_string(__name__, template_path).decode("UTF-8")

    file_name = "{}/{}/time_k_best_t{}.sh".format(replacement['linnea_jobscripts_path'], replacement["name"], replacement["threads"])
    with open(file_name, "wt", encoding='utf-8') as output_file:
        print("Writing", file_name)
        output_file.write(template_str.format(runner_name="runner_k_best_t{}".format(replacement["threads"]), output_subdir="k_best", **replacement))

def generate_code_scripts(replacement):

    template_path = "jobscripts/templates/generate_code.sh"
    template_str = pkg_resources.resource_string(__name__, template_path).decode("UTF-8")

    for args, ref in [("", ""), ("-f", "_ref")]:
        replacement_copy = replacement.copy()

        replacement_copy["args"] = args
        replacement_copy["ref"] = ref

        file_name = "{}/{}/generate_code{}.sh".format(replacement_copy['linnea_jobscripts_path'],
                                                       replacement_copy["name"], ref)
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
        "flag_node":            "-m ",
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
        "flag_model":           "-C",
        "flag_exclusive":       "--exclusive",
        "flag_node":            "--nodelist=",
        "var_array_idx":        "SLURM_ARRAY_TASK_ID",
        "string_array_idx":     "%a"
        }
    }

def generate_scripts(experiment, number_of_experiments,  num_threads):

    experiment_configuration = config.experiment_configuration

    scheduler = experiment_configuration["scheduler"]

    if experiment_configuration["time"]["exclusive"]:
        experiment_configuration["time"]["spec_exclusive"] = "#{directive} {flag_exclusive}".format(**scheduler_vars[scheduler])
    else:
        experiment_configuration["time"]["spec_exclusive"] = ""

    if experiment_configuration["time"]["node"]:
        experiment_configuration["time"]["spec_node"] = "#{directive} {flag_node}".format(**scheduler_vars[scheduler]) + experiment_configuration["time"]["node"]
    else:
        experiment_configuration["time"]["spec_node"] = ""

    for mode in ["time", "generate"]:
        if scheduler == "LSF":
            experiment_configuration[mode]["lsf_arrayjob"] = "[1-{}]".format(number_of_experiments)
            experiment_configuration[mode]["slurm_arrayjob"] = ""
        elif scheduler == "SLURM":
            experiment_configuration[mode]["lsf_arrayjob"] = ""
            experiment_configuration[mode]["slurm_arrayjob"] = "#SBATCH --array=1-{}".format(number_of_experiments)
        else:
            pass # TODO throw exception?

        experiment_configuration[mode]['name'] = experiment


    dirname = "{}/{}/".format(experiment_configuration['path']['linnea_jobscripts_path'], experiment)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    time_configuration = {**experiment_configuration['time'], **experiment_configuration['path'], **experiment_configuration['version'], **scheduler_vars[scheduler]}
    generate_configuration = {**experiment_configuration['generate'], **experiment_configuration['path'], **experiment_configuration['version'], **scheduler_vars[scheduler]}

    for threads in num_threads:
        time_configuration_copy = time_configuration.copy()
        time_configuration_copy["threads"] = threads
        time_execution_scripts(time_configuration_copy)
        time_k_best_script(time_configuration_copy)

    time_generation_script(time_configuration)
    generate_code_scripts(generate_configuration)

if __name__ == '__main__':
    generate_scripts()
