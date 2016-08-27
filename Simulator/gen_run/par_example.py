import lib_xml as lxml
import os 
import global_par as g
import sys

#---------------------------------------------------------
# define the basic directories in global_par.py

dir_DB=g.dir_proj+'DB/'
g_runID='20160818_140GHz_TQU_7'

#---------------------------------------------------------
# set the batch run mode
debug='regular'

# <log>
runID=g_runID
dir_simulator=g.dir_simulator
machine='kekcc'
silent='Y'

# <basicpar>
nside=256
pixelmapio='N'
run_type='sim'
coord='E2G'

# new addition
#TQU='T_bandpassmismatch'
TQU='TQU'
gen_tod='N'
# </basicpar>

# <simulations>
file_input_maps='/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM_new/processed_maps/140GHz/thermo/fits/cmb_map_band40_b32arcmin_ud512_thermo'
file_input_maps2='/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM_new/processed_maps/280GHz/thermo/fits/cmb_map_band30_FG_sidelobe_50dB_b1200arcmin_ud512_thermo'

# new addition
file_input_map_dust='/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM_new/processed_maps/280GHz/thermo/fits/cmb_map_band30_FG_sidelobe_50dB_b1200arcmin_ud512_thermo'
file_input_map_synch='/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM_new/processed_maps/280GHz/thermo/fits/cmb_map_band30_FG_sidelobe_50dB_b1200arcmin_ud512_thermo'

file_input_noise=dir_DB+'sim_nonoise.txt'
file_input_fpdb=dir_DB+'/fp_db/LB_HFW_example_20160813.db'
file_input_muellermatrix=dir_DB+'pb1_MuellerMatrix_default.db'

# new addition
file_input_bandpassmismatch=dir_DB+'/sys_db/test_alphaD.txt'
# </simulations>

# <database>
#file_fpdb_mmin=dir_DB+'/fp_db/polang_random/LB_HFW_example_sigma1degs.db'
file_fpdb_mmin=dir_DB+'/fp_db/LB_HFW_example_20160813.db'
file_muellermatrix=dir_DB+'pb1_MuellerMatrix_default.db'
#file_relgain=dir_DB+'pb1_simgain_default.db'
#file_boloid=dir_DB+'/fp_db/LB_HFW_example_v2.db'
# </database>

# <filtering_choice>
filter_choice='poly'
poly=-1
#file_noisefft=dir_DB+'pb1_observation_fake_noHWP.db'
# </filtering_choice> 

# <output_file>
dir_simedmap=g.dir_proj+'/SimedMaps/RunLog/'+runID
dir_combinedmap=dir_simedmap
Tn_map='map_Tn.fits'
Td_map='map_Td.fits'
AA_map='map_AA.fits'
BB_map='map_BB.fits'
AB_map='map_AB.fits'
Ad_map='map_Ad.fits'
Bd_map='map_Bd.fits'
# </output_file>

# <data_selection>
dir_MMrunID=g.dir_proj+'/SimedMaps/RunLog/'+runID
#file_flag_pixel=dir_DB+'/beamprm_20120530_031419_hwp112.5.db_simin_flag.db'
db_gain=dir_DB+'/gain_db/relgain_flat_bias_0.db'
gain_type='ideal'
db_ces='/group/cmb/litebird/simdata/Scans/dataout/LB_L2_20131226_samplerate_SCANSPEC_1440min_65.0degs_93.0min_30.0degs_0.1rpm_365day_nside256_10Hz/ptg/LBPTG_latlon_eclip.db'
db_ces2='/home/cmb/tmatsumu/develop/LiteBIRD/projects/LB_SYSPL_v3.1/pyScans/dataout/LB_L2_20140205_samplerate_SCANSPEC_1440min_65.0degs_93.0min_55.0degs_0.1rpm_365day_nside256_10Hz/ptg/LBPTG_latlon_eclip.db'
# </data_selection>


#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
# prepare for the directories
if os.path.exists(dir_MMrunID):
    print '[par.py] ', 'runID directory exists' 
if not os.path.exists(dir_MMrunID):
    print '[par.py] ', 'runID directory does not exist as'+dir_MMrunID
    print '[par.py] ', 'creating '+dir_MMrunID 
    os.system('mkdir -p '+dir_MMrunID)

