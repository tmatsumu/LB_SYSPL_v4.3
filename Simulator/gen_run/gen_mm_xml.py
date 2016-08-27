import sys
import lib_xml as lxml
#import lib_mapmaker as lmm
import importlib
param_module=sys.argv[1]
params = importlib.import_module(param_module)

fname_xml = sys.argv[2]
lxml.write_xml(fname_xml, params)

sys.exit()
