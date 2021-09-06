import enum
import importlib
import os.path
import json
import math
import logging

logger = logging.getLogger('config')

_LOCAL_CONFIG_FILE = 'linnea_config.json'
_GLOBAL_CONFIG_FILE = os.path.expandvars('$HOME/linnea_config.json')

class LanguageNotSet(Exception):
    pass

class OutputPathNotSet(Exception):
    pass

class DataTypeNotSet(Exception):
    pass

class UnsupportedDataType(Exception):
    pass

class LanguageOptionNotImplemented(Exception):
    pass

class DirectoryDoesNotExist(Exception):
    pass

class CDataType(enum.Enum):
    float = 0
    double = 1

class JuliaDataType(enum.Enum):
    Float32 = 0
    Float64 = 1

class Language(enum.Enum):
    C = 0
    Julia = 1
    Matlab = 2
    Cpp = 3

class CppLibrary(enum.Enum):
    Blaze = 0
    Eigen = 1
    Armadillo = 2

class GraphStyle(enum.Enum):
    full = 0
    simple = 1
    minimal = 2

# defining variables

c = False
julia = False
matlab = False

float = False
double = False

float32 = False
float64 = False

data_type_string = None
blas_data_type_prefix = None
filename_extension = None
comment = None

language = None
output_name = None
merging_branches = False
dead_ends = False
pruning_factor = -1
algorithms_limit = -1
generate_graph = False
generate_derivation = False
generate_code = False
generate_experiments = False
verbosity = -1
time_limit = -1
graph_style = None
experiment_configuration = dict()
output_code_path = None
results_path = None

def set_language(_language):
    global language, filename_extension, comment, julia, c, matlab
    language = _language
    if _language is Language.C:
        c = True
        julia = False
        matlab = False
        filename_extension = ".c"
        comment = "// "
    elif _language is Language.Julia:
        julia = True
        c = False
        matlab = False
        filename_extension = ".jl"
        comment = "# "
    else:
        raise LanguageOptionNotImplemented()

def import_collections():
    global language
    if language:
        return importlib.import_module("linnea.kernels.{}.collections".format(language.name.lower()))
    else:
        raise LanguageNotSet()

def set_data_type(data_type):
    global data_type_string
    data_type_string = data_type.name
    global language
    if language is Language.C:
        global blas_data_type_prefix, float, double
        if data_type is CDataType.float:
            blas_data_type_prefix = "s"
            float = True
            double = False
        elif data_type is CDataType.double:
            blas_data_type_prefix = "d"
            double = True
            float = False
        else:
            raise UnsupportedDataType()
    elif language is Language.Julia:
        global float32, float64
        if data_type is JuliaDataType.Float32:
            float32 = True
            float64 = False
        elif data_type is JuliaDataType.Float64:
            float64 = True
            float32 = False
        else:
            raise UnsupportedDataType()

def set_merging_branches(merging):
    global merging_branches
    merging_branches = merging

def set_dead_ends(de):
    global dead_ends
    dead_ends = de

def set_pruning_factor(factor):
    global pruning_factor
    pruning_factor = factor

def set_algorithms_limit(n):
    global algorithms_limit
    algorithms_limit = n

def set_generate_graph(generate):
    global generate_graph
    generate_graph = generate

def set_output_name(name):
    global output_name
    output_name = name

def set_generate_derivation(generate):
    global generate_derivation
    generate_derivation = generate

def set_generate_code(generate):
    global generate_code
    generate_code = generate

def set_generate_experiments(generate):
    global generate_experiments
    generate_experiments = generate

def set_verbosity(level):
    global verbosity
    verbosity = level

def set_time_limit(limit):
    global time_limit
    if isinstance(limit, str):
        if limit == "inf":
            time_limit = math.inf
        else:
            time_limit = float(limit)
    else:    
        time_limit = limit

def set_graph_style(style):
    global graph_style
    graph_style = style

def init():
    pass
    # TODO This function used to do something, but right now it's not needed anymore. Remove?

def clear_all():
    from linnea import temporaries
    from linnea.algorithm_generation import partitioning
    from linnea.kernels.julia import collections
    temporaries.clear()
    partitioning.clear()
    collections.clear()

def create_path(path):

    path = os.path.abspath(os.path.expanduser(os.path.expandvars(path)))
    if not os.path.exists(path):
        os.makedirs(path)
    return path

def check_path(path):

    path = os.path.abspath(os.path.expanduser(os.path.expandvars(path)))
    if not os.path.exists(path):
        raise DirectoryDoesNotExist(path)
    else:
        return path

def set_output_code_path(path):
    global output_code_path
    output_code_path = create_path(path)

def set_results_path(path):
    global results_path
    results_path = create_path(path)

# setting default values

set_language(Language.Julia)
set_data_type(JuliaDataType.Float64)
set_merging_branches(True)
set_dead_ends(True)
set_pruning_factor(1.0)
set_algorithms_limit(100)
set_generate_graph(False)
set_generate_derivation(False)
set_generate_code(True)
set_generate_experiments(False)
set_output_code_path('.')
set_verbosity(1)
set_time_limit(60)
set_graph_style(GraphStyle.full)

def load_config(config_file=''):

    global experiment_configuration

    if config_file:
        config_file = os.path.expandvars(config_file)
        if not os.path.exists(config_file):
            msg = "Specified config file '{}' not found.".format(config_file)
            raise FileNotFoundError(msg)
    else:
        if os.path.exists(_LOCAL_CONFIG_FILE):
            config_file= _LOCAL_CONFIG_FILE
        elif os.path.exists(_GLOBAL_CONFIG_FILE):
            config_file = _GLOBAL_CONFIG_FILE

    if config_file:
        logger.info('Found config file: {}'.format(os.path.abspath(config_file)))
        with open(config_file) as jsonfile:
            settings = globals()

            configuration = json.load(jsonfile)

            # set up configuration for normal use
            for key, value in configuration['main'].items():
                if key == 'language':
                    set_language(Language[value])
                elif key == 'c_data_type':
                    set_data_type(CDataType[value])
                elif key == 'julia_data_type':
                    set_data_type(JuliaDataType[value])
                elif key == 'output_code_path':
                    set_output_code_path(value)
                elif key == 'time_limit':
                    set_time_limit(value)
                elif key == 'graph_style':
                    set_graph_style(GraphStyle[value])
                elif key in settings and not key.startswith('_'):
                    settings[key] = value
                else:
                    error_msg = 'Unknown setting: {}'.format(key)
                    raise KeyError(error_msg)

            if 'experiments' in configuration.keys():

                # create a configuration dictionaries for experiments
                for key, value in configuration['experiments']['path'].items():
                    if key in ['linnea_src_path', 'linnea_lib_path', 'linnea_julia_path', 'linnea_virtualenv_path']:
                        configuration['experiments']['path'][key] = check_path(value)
                    elif key in ['linnea_output_path', 'linnea_results_path', 'linnea_jobscripts_path']:
                        configuration['experiments']['path'][key] = create_path(value)
                        if key == 'linnea_results_path':
                            set_results_path(value)
                    else:
                        error_msg = 'Unknown setting: {}'.format(key)
                        raise KeyError(error_msg)

                configuration['experiments']['path']['output_code_path'] = output_code_path

                experiment_configuration = configuration['experiments']
    else:
        logger.info('No config file found. Using default settings.')


load_config()