#!/bin/sh

arcsec=$1

dir_viewer=/global/homes/t/tmatsumu/develop/PBI/repo_test/PB1_NTP_develop/src_viewer/

dir_db=/scratch/scratchdirs/tmatsumu/sim/PB1_NTP/DB/
dbname=beamprm_20120530_031419_hwp112.5.db_simin_ave_fact1.db

#python gen_FPDB_simplifieddiffptg.py $dir_db $dbname 5.
python gen_FPDB_simplifieddiffptg.py $dir_db $dbname beamprm_20120530_031419_hwp112.5.db_simin_${arcsec}arcsec_simplified_diffptg.db $arcsec
cp ${dir_db}/beamprm_20120530_031419_hwp112.5.db_simin_flag_fact1.db ${dir_db}beamprm_20120530_031419_hwp112.5.db_simin_flag_${arcsec}arcsec_simplified_diffptg.db
cp ${dir_db}/beamprm_20120530_031419_hwp112.5.db_simin_ave_fact1.db ${dir_db}beamprm_20120530_031419_hwp112.5.db_simin_ave_${arcsec}arcsec_simplified_diffptg.db
exit

python ${dir_viewer}viewFPDB_polangvec.py $dir_db $beamprm_20120530_031419_hwp112.5.db_simin_${arcsec}arcsec_simplified_diffptg.db beamprm_20120530_031419_hwp112.5.db_simin_flag_${arcsec}arcsec_simplified_diffptg.db 5e-4
mv ${dir_db}${dbname}_simplified_diffptg.db_pixel.png ~/www/tmp
mv ${dir_db}${dbname}_simplified_diffptg.db_pixeloverdiffptg.png ~/www/tmp
mv ${dir_db}${dbname}_simplified_diffptg.db_diffptg.png ~/www/tmp
chmod 777 ~/www/tmp/*
