import pkg_resources
import json
import os.path
import itertools

# this causes a problem because it tries to read the config.json in the current directory.
# from linnea.config import Strategy

def main():

    replacement = dict()
    
    config_file = "config.json"
    if os.path.exists(config_file):
        with open(config_file) as jsonfile:
            for key, value in json.load(jsonfile).items():
                if key == "exclusive":
                    if value:
                        replacement["exclusive"] = "#BSUB -x                     # exclusive access"
                    else:
                        replacement["exclusive"] = ""
                else:
                    replacement[key] = value    

    template_path = "templates/time_generation.sh"
    template_str = pkg_resources.resource_string(__name__, template_path).decode("UTF-8")

    file_name = "time_generation_test.sh"
    for strategy, merging in itertools.product(["c", "e"], [True, False]):
        replacement_copy = replacement.copy()

        file_name_parts = ["time_gen"]

        file_name_parts.append(strategy)
        replacement_copy["strategy"] = strategy
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

        file_name_parts.append(".sh")
        file_name = "_".join(file_name_parts)
        with open(file_name, "wt", encoding='utf-8') as output_file:
            print("Writing", file_name)
            output_file.write(template_str.format(**replacement_copy))


if __name__ == '__main__':
    main()