from numpy import *
import sys
import glob
import ReadMapMakeXml as rxml

def Gen_scansetdirs(xml_filename):
    
    # xml_filename = sys.argv[1]
    input_xml = rxml.Get_Mapmake_Inputs(xml_filename)

    fname_ptg = input_xml["input_ptg"]
        
#    date_i = input_xml["date_i"]
#    date_f = input_xml["date_f"]
#    fsample = input_xml["fsample"]
#    nside = input_xml["nside"]
#    run_type = input_xml["run_type"]
#    fname_map = input_xml["input_maps"]
#    fname_noise = input_xml["input_noise"]
#    fname_fpdb = input_xml["fpdb_file"]
#    fname_detectors = input_xml["detectors"]
#    fname_HWPangles = input_xml["HWPangles"]
#    fname_scanset = input_xml["scan_set"]
#    poly = input_xml["poly"]
#    dir_simedmap = input_xml["dir_simedmap"]
    
#    fnameout_Tn = input_xml["Tn"]
#    fnameout_Tn = input_xml["Td"]
#    fnameout_Tn = input_xml["AA"]
#    fnameout_Tn = input_xml["BB"]
#    fnameout_Tn = input_xml["AB"]
#    fnameout_Tn = input_xml["Ad"]
#    fnameout_Tn = input_xml["Bd"]
    
    fileNames = glob.glob(fname_ptg+"/observation_*")
#    print fname_ptg
    nb = len(fileNames)
        
    fileNames_date = zeros(nb,int)
    for i in range(0,nb): fileNames_date[i] = int(fileNames[i][89:]) 
    
    ind = where((fileNames_date >= date_i) & (fileNames_date <= date_f))
    
    nb_obs = len(ind[0])
    dir_obs = []
    for i in range(0,nb_obs): dir_obs.append(fname_ptg+"/observation_"+str(fileNames_date[ind[0][i]]))
    
    file_input_pointing = []
    for i in range(0,nb_obs):
        fileNames = glob.glob(dir_obs[i]+'/scan*')
        nb = len(fileNames)
        #    print fileNames
        for j in range(0,nb): 
            file_input_pointing.append(fileNames[j])
            
    nb_tmp = nb_obs*nb
    for i in range(0,nb_tmp): print file_input_pointing[i]

    input_xml["file_input_pointing"] = file_input_pointing

    return input_xml
