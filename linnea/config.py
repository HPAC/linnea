import enum
import importlib
import os.path
import json
import math

_LOCAL_CONFIG_FILE = 'config.json'
_GLOBAL_CONFIG_FILE = os.path.expandvars('$HOME/.linnea.config.json')

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

class Strategy(enum.Enum):
    constructive = 0
    exhaustive = 1

class GraphStyle(enum.Enum):
    full = 0
    simple = 1
    minimal = 2

# defining variables

language = None
output_path = None
output_name = None
merging_branches = False
dead_ends = False
algorithms_limit = -1
generate_graph = False
generate_derivation = False
generate_code = False
generate_experiments = False
strategy = None
verbosity = -1
solution_nodes_limit = -1
iteration_limit = -1
graph_style = None
experiment_configuration = dict()
linnea_install_path = '$HOME/Linnea'
linnea_src_path = '$HOME/Linnea/src/linnea'
linnea_lib_path = '$HOME/Linnea/lib'
linnea_output_path = '$HOME/Linnea/output'
linnea_results_path = '$HOME/Linnea/output/results'
linnea_jobscripts_path = '$HOME/Linnea/output/jobscripts'
linnea_code_path = '$HOME/Linnea/output/code'
linnea_virtualenv_path = '$HOME/Linnea/linnea.venv'
linnea_julia_path = '$HOME/Linnea/src/julia'
linnea_job_log_path = '$HOME/Linnea/output/job_logs'

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

def set_strategy(_strategy):
    global strategy
    strategy = _strategy

def set_merging_branches(merging):
    global merging_branches
    merging_branches = merging

def set_dead_ends(de):
    global dead_ends
    dead_ends = de

def set_algorithms_limit(n):
    global algorithms_limit
    algorithms_limit = n

def set_generate_graph(generate):
    global generate_graph
    generate_graph = generate

def set_output_name(name):
    global output_name
    output_name = name

def set_output_path(path):
    global output_path
    output_path = os.path.abspath(os.path.expanduser(os.path.expandvars(path)))
    if not os.path.exists(output_path):
        os.makedirs(output_path)

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

def set_solution_nodes_limit(limit):
    global solution_nodes_limit
    if isinstance(limit, str):
        if limit == "inf":
            solution_nodes_limit = math.inf
        else:
            solution_nodes_limit = int(limit)
    else:    
        solution_nodes_limit = limit

def set_iteration_limit(limit):
    global iteration_limit
    if isinstance(limit, str):
        if limit == "inf":
            iteration_limit = math.inf
        else:
            iteration_limit = int(limit)
    else:    
        iteration_limit = limit

def set_graph_style(style):
    global graph_style
    graph_style = style

def init():
    pass
    # TODO This function used to do something, but right now it's not needed anymore. Remove?

def clear_all():
    from linnea import temporaries
    from linnea.derivation import partitioning
    temporaries.clear()
    partitioning.clear()


# setting default values

set_language(Language.Julia)
set_data_type(JuliaDataType.Float64)
set_merging_branches(True)
set_dead_ends(True)
set_algorithms_limit(100)
set_generate_graph(False)
set_generate_derivation(False)
set_generate_code(True)
set_generate_experiments(False)
set_strategy(Strategy.constructive)
set_output_path('.')
set_verbosity(1)
set_solution_nodes_limit(math.inf)
set_iteration_limit(100)
set_graph_style(GraphStyle.full)


def load_config():

    global experiment_configuration
    config_file= ''
    if os.path.exists(_LOCAL_CONFIG_FILE):
        config_file= _LOCAL_CONFIG_FILE
    elif os.path.exists(_GLOBAL_CONFIG_FILE):
        config_file = _GLOBAL_CONFIG_FILE

    if config_file:
        with open(config_file) as jsonfile:
            settings = globals()

            configuration = json.load(jsonfile)

            # make sure all paths are properly expanded
            configuration['path'] = {k: os.path.expandvars(v) for k, v in configuration['path'].items()}

            # create a configuration dictionary for normal use and the experiments
            configuration['main'] = {**configuration['main'], **configuration['path']}
            configuration['experiments']['time'] = {**configuration['experiments']['time'], **configuration['path']}
            configuration['experiments']['generate'] = {**configuration['experiments']['generate'], **configuration['path']}

            # set up configuration for normal use
            for key, value in configuration['main'].items():
                if key == 'language':
                    set_language(Language[value])
                elif key == 'c_data_type':
                    set_data_type(CDataType[value])
                elif key == 'julia_data_type':
                    set_data_type(JuliaDataType[value])
                elif key == 'strategy':
                    set_strategy(Strategy[value])
                elif key == 'linnea_code_path':
                    set_output_path(value)
                elif key == 'solution_nodes_limit':
                    set_solution_nodes_limit(value)
                elif key == 'iteration_limit':
                    set_iteration_limit(value)
                elif key == 'graph_style':
                    set_graph_style(GraphStyle[value])
                elif key in settings and not key.startswith('_'):
                    settings[key] = value
                else:
                    error_msg = 'Unknown setting: {}'.format(key)
                    raise KeyError(error_msg)

            # set up configuration for experiments
            if 'exclusive' in configuration['experiments']['time'].keys():
                if configuration['experiments']['time']['exclusive']:
                    configuration['experiments']['time']['exclusive'] = '#BSUB -x                     # exclusive access'
                else:
                    configuration['experiments']['time']['exclusive'] = ''

            experiment_configuration = configuration['experiments']

load_config()