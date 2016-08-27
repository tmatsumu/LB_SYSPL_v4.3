import matsumulib as mylib 

#dir_wrk='/raid/users/tmatsumu/work/develop/LiteBIRD/ScanStrategy/pyScans/'
#dir_dataout='/raid/users/tmatsumu/data_sim/LiteBIRD/Scans/'
dir_wrk='/home/cmb/tmatsumu/develop/LiteBIRD/projects/20130924_LBBasicMM/pyScans/'
dir_dataout='/group/cmb/litebird/simdata/Scans/'

path_inputparams=dir_wrk+'input_params/'
path_src=dir_wrk+'src/'
path_out=dir_dataout+'dataout/'

today_julian = mylib.convert_Gregorian2Julian( 2020, 3, 20, 3, 50, 0)
option_gen_ptg = True
