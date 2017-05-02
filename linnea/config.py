import enum
import importlib
import os.path

class LanguageNotSet(Exception):
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

language = None

output_path = os.path.expanduser("~/Promotion/code/workspace")
if not os.path.exists(output_path):
    raise DirectoryDoesNotExist(output_path)

experiments_path = os.path.expanduser("~/Promotion/code/experiments")
if not os.path.exists(experiments_path):
    raise DirectoryDoesNotExist(experiments_path)

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

def set_language(_language):
    global language, filename_extension, comment
    language = _language
    if _language is Language.C:
        global production_code_fragment_type, c
        c = True
        filename_extension = ".c"
        comment = "// "
    elif _language is Language.Julia:
        global julia
        julia = True
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
        global blas_data_type_prefix
        if data_type is CDataType.float:
            blas_data_type_prefix = "s"
            global float
            float = True
        elif data_type is CDataType.double:
            blas_data_type_prefix = "d"
            global double
            double = True
        else:
            raise UnsupportedDataType()
    elif language is Language.Julia:
        if data_type is JuliaDataType.Float32:
            global float32
            float32 = True
        elif data_type is JuliaDataType.Float64:
            global float64
            float64 = True
        else:
            raise UnsupportedDataType()        

def init():
    from .algebra import property_DNs
    property_DNs._init()
