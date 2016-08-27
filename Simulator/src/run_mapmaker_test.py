import sys
import lib_mapmaker as lmm
import numpy as np

def npojt2dict(fname):
    outputs = {}
    inputs = np.load(fname)
    outputs['SimInput']=inputs[()]["SimInput"]
    outputs['fpdb_list_simgen']=inputs[()]["fpdb_list_simgen"]
    outputs['fpdb_list_mmin']=inputs[()]["fpdb_list_mmin"]
    outputs['NoiseInput']=inputs[()]["NoiseInput"]
    outputs['pix_list']=inputs[()]["pix_list"]
#    outputs['SysInputs']=inputs[()]["SysInputs"]
    outputs['poly']=inputs[()]["poly"]
    outputs['nside']=inputs[()]["nside"]
    outputs['runtime_init']=inputs[()]["runtime_init"]
    outputs['out_dir']=inputs[()]["out_dir"]
    outputs['out_dir_ptg']=inputs[()]["out_dir_ptg"] 
    outputs['ptg_package']=inputs[()]["ptg_package"]
    outputs['ptg_idx']=inputs[()]["ptg_idx"]
    outputs['tmp']=inputs[()]["tmp"]
    outputs['muellermatrix']=inputs[()]["muellermatrix"]
    outputs['relgain']=inputs[()]["relgain"]
    del(inputs)
    return(outputs)

fname_inputs = sys.argv[1]
input = npojt2dict(fname_inputs)
#lmm.main_mapmaker(input)

pix_list = input["pix_list"]
print pix_list
