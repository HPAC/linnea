import sys
import os.path

_test_path = os.path.dirname(os.path.abspath(__file__))
_linnea_path = os.path.abspath(os.path.join(_test_path, '..'))
sys.path.insert(0, _linnea_path)
