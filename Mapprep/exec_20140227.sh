#!/bin/sh

# filename_map_in filename_map_out nside_out beam [arcmins] supression

nside_out=512
dir_h=/home/cmb/tmatsumu/develop/LiteBIRD/projects/LB_SYSPL_v3.2/Mapprep/

band_name=$1
unit=$2
beam_arcmin_sidelobe=1200

if [[ ${band_name}GHz = "60GHz" ]];then
    echo ${band_name}GHz
    beam_arcmin=75
    band=band10
    dir=/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM_new/000uK/lensed/r00001/k10/0051-0100/PSM_OUTPUT/skyinbands/PSM_IDEAL/${band}/
fi

if [[ ${band_name}GHz = "78GHz" ]];then
    echo ${band_name}GHz
    beam_arcmin=60
    band=band28
    dir=/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM_new/000uK/lensed/r00001/k10/0051-0100/PSM_OUTPUT/skyinbands/PSM_IDEAL/${band}/
fi

if [[ ${band_name}GHz = "100GHz" ]];then
    echo ${band_name}GHz
    beam_arcmin=45
    band=band50
    dir=/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM_new/000uK/lensed/r00001/k10/0051-0100/PSM_OUTPUT/skyinbands/PSM_IDEAL/${band}/
fi

if [[ ${band_name}GHz = "140GHz" ]];then
    echo ${band_name}GHz
    beam_arcmin=32
    band=band40
    dir=/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM_new/000uK/lensed/r00001/k10/0101-0150/PSM_OUTPUT/skyinbands/PSM_IDEAL/${band}/
fi

if [[ ${band_name}GHz = "195GHz" ]];then
    echo ${band_name}GHz
    beam_arcmin=23
    band=band45
    dir=/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM_new/000uK/lensed/r00001/k10/0151-0200/PSM_OUTPUT/skyinbands/PSM_IDEAL/${band}/
fi

if [[ ${band_name}GHz = "280GHz" ]];then
    echo ${band_name}GHz
    beam_arcmin=16
    band=band30
    dir=/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM_new/000uK/lensed/r00001/k10/0251-0300/PSM_OUTPUT/skyinbands/PSM_IDEAL/${band}/
fi

dir_out=/group/cmb/litebird/simdata/Maps/20131001_linux_711_LSM_new/processed_maps/${band_name}GHz/${unit}
mkdir -p $dir_out

#+++++++++++++++++++++++
# note:
# gen_sidelobemap_dipole inputmappath(no .fits) outputFilename(no .fits) nside_out band_GHz beam_arcmin dB option_dipole(y/n) unit_option(thermo/RJ) option_FG(y/n) file_synch file_dust

#+++++++++++++++++++++++
# raw map

# CMB
filename_in=cmb_map_$band
filename_out=${filename_in}_b0arcmin_ud1024
bsub -o $dir_h/log/log1.txt -q e python gen_sidelobemap_dipole.py $dir/$filename_in $dir_out/$filename_out 1024 ${band_name} 0. 0 n $unit n

# Synch
filename_in=synchrotron_map_$band
filename_out=${filename_in}_b0arcmin_ud1024
bsub -o $dir_h/log/log2.txt -q e python gen_sidelobemap_dipole.py $dir/$filename_in $dir_out/$filename_out 1024 ${band_name} 0. 0 n $unit n

# Dust
filename_in=thermaldust_map_$band
filename_out=${filename_in}_b0arcmin_ud1024
bsub -o $dir_h/log/log3.txt -q e python gen_sidelobemap_dipole.py $dir/$filename_in $dir_out/$filename_out 1024 ${band_name} 0. 0 n $unit n

#+++++++++++++++++++++++
# main beam map, no dipole

# CMB
filename_in=cmb_map_$band
filename_out=${filename_in}_b${beam_arcmin}arcmin_ud${nside_out}
bsub -o $dir_h/log/log4.txt -q e python gen_sidelobemap_dipole.py $dir/$filename_in $dir_out/$filename_out $nside_out ${band_name} $beam_arcmin 0 n $unit n

# Synch
filename_in=synchrotron_map_$band
filename_out=${filename_in}_b${beam_arcmin}arcmin_ud${nside_out}
bsub -o $dir_h/log/log5.txt -q e python gen_sidelobemap_dipole.py $dir/$filename_in $dir_out/$filename_out $nside_out ${band_name} $beam_arcmin 0 n $unit n

# Dust
filename_in=thermaldust_map_$band
filename_out=${filename_in}_b${beam_arcmin}arcmin_ud${nside_out}
bsub -o $dir_h/log/log6.txt -q e python gen_sidelobemap_dipole.py $dir/$filename_in $dir_out/$filename_out $nside_out ${band_name} $beam_arcmin 0 n $unit n
#++++++++++++++++++++++++
# main beam map, w/ Dipole

# CMB+Dipole
filename_in=cmb_map_$band
filename_out=${filename_in}_Dipole_b${beam_arcmin}arcmin_ud${nside_out}
bsub -o $dir_h/log/log7.txt -q e python gen_sidelobemap_dipole.py $dir/$filename_in $dir_out/$filename_out $nside_out ${band_name} $beam_arcmin 0 y $unit n

# Synch+Dipole
filename_in=synchrotron_map_$band
filename_out=${filename_in}_Dipole_b${beam_arcmin}arcmin_ud${nside_out}
bsub -o $dir_h/log/log8.txt -q e python gen_sidelobemap_dipole.py $dir/$filename_in $dir_out/$filename_out $nside_out ${band_name} $beam_arcmin 0 y $unit n

# Dust+Dipole
filename_in=thermaldust_map_$band
filename_out=${filename_in}_Dipole_b${beam_arcmin}arcmin_ud${nside_out}
bsub -o $dir_h/log/log9.txt -q e python gen_sidelobemap_dipole.py $dir/$filename_in $dir_out/$filename_out $nside_out ${band_name} $beam_arcmin 0 y $unit n
#++++++++++++++++++++++++

# CMB+Synch+Dust
filename_cmb=cmb_map_$band
filename_dust=thermaldust_map_$band
filename_synch=synchrotron_map_$band
filename_out=${filename_cmb}_FG_b${beam_arcmin}arcmin_ud${nside_out}
bsub -o $dir_h/log/log10.txt -q e python gen_sidelobemap_dipole.py $dir/$filename_cmb $dir_out/$filename_out $nside_out ${band_name} $beam_arcmin 0 n $unit y $dir/$filename_synch $dir/$filename_dust

# CMB+Synch+Dust+Dipole
filename_cmb=cmb_map_$band
filename_dust=thermaldust_map_$band
filename_synch=synchrotron_map_$band
filename_out=${filename_cmb}_FG_Dipole_b${beam_arcmin}arcmin_ud${nside_out}
bsub -o $dir_h/log/log11.txt -q e python gen_sidelobemap_dipole.py $dir/$filename_cmb $dir_out/$filename_out $nside_out ${band_name} $beam_arcmin 0 y $unit y $dir/$filename_synch $dir/$filename_dust

#+++++++++++++++++++++++
#+++++++++++++++++++++++
# sidelobe map

# Synch+Dipole, 0-100dB
filename_in=synchrotron_map_${band}
for((i=0 ; i<=100 ; i+=10))
do
    filename_out=${filename_in}_Dipole_sidelobe_${i}dB_b${beam_arcmin_sidelobe}arcmin_ud${nside_out}
    bsub -o $dir_h/log/log12_${i}.txt -q e python gen_sidelobemap_dipole.py $dir/$filename_in $dir_out/$filename_out $nside_out ${band_name} $beam_arcmin_sidelobe -$i y $unit n
done

#+++++++++++++++++++++++

# Dust+Dipole, 0-100dB
filename_in=thermaldust_map_${band}
for((i=0 ; i<=100 ; i+=10))
do 
    filename_out=${filename_in}_Dipole_sidelobe_${i}dB_b${beam_arcmin_sidelobe}arcmin_ud${nside_out}
    bsub -o $dir_h/log/log13_${i}.txt -q e python gen_sidelobemap_dipole.py $dir/$filename_in $dir_out/$filename_out $nside_out ${band_name} $beam_arcmin_sidelobe -$i y $unit n
done

#+++++++++++++++++++++++

# CMB+Dipole, 0-100dB 
filename_in=cmb_map_$band
for((i=0 ; i<=100 ; i+=10))
do 
    filename_out=${filename_in}_Dipole_sidelobe_${i}dB_b${beam_arcmin_sidelobe}arcmin_ud${nside_out}
    bsub -o $dir_h/log/log14_${i}.txt -q e python gen_sidelobemap_dipole.py $dir/$filename_in $dir_out/$filename_out $nside_out ${band_name} $beam_arcmin_sidelobe -$i y $unit n
done

#+++++++++++++++++++++++
#+++++++++++++++++++++++

# Synch, 0-100dB
filename_in=synchrotron_map_${band}
for((i=0 ; i<=100 ; i+=10))
do 
    filename_out=${filename_in}_sidelobe_${i}dB_b${beam_arcmin_sidelobe}arcmin_ud${nside_out}
    bsub -o $dir_h/log/log15_${i}.txt -q e python gen_sidelobemap_dipole.py $dir/$filename_in $dir_out/$filename_out $nside_out ${band_name} $beam_arcmin_sidelobe -$i n $unit n
done

#+++++++++++++++++++++++

# Dust, 0-100dB
filename_in=thermaldust_map_${band}
for((i=0 ; i<=100 ; i+=10))
do 
    filename_out=${filename_in}_sidelobe_${i}dB_b${beam_arcmin_sidelobe}arcmin_ud${nside_out}
    bsub -o $dir_h/log/log16_${i}.txt -q e python gen_sidelobemap_dipole.py $dir/$filename_in $dir_out/$filename_out $nside_out ${band_name} $beam_arcmin_sidelobe -$i n $unit n
done

#+++++++++++++++++++++++

# CMB, 0-100dB
filename_in=cmb_map_$band
for((i=0 ; i<=100 ; i+=10))
do 
    filename_out=${filename_in}_sidelobe_${i}dB_b${beam_arcmin_sidelobe}arcmin_ud${nside_out}
    bsub -o $dir_h/log/log17_${i}.txt -q e python gen_sidelobemap_dipole.py $dir/$filename_in $dir_out/$filename_out $nside_out ${band_name} $beam_arcmin_sidelobe -$i n $unit n
done

#+++++++++++++++++++++++
#+++++++++++++++++++++++

# CMB+Synch+Dust
filename_cmb=cmb_map_$band
filename_dust=thermaldust_map_$band
filename_synch=synchrotron_map_$band
for((i=0 ; i<=100 ; i+=10))
do 
    filename_out=${filename_cmb}_FG_sidelobe_${i}dB_b${beam_arcmin_sidelobe}arcmin_ud${nside_out}
    bsub -o $dir_h/log/log18_${i}_${band_name}_${unit}.txt -q e python gen_sidelobemap_dipole.py $dir/$filename_cmb $dir_out/$filename_out $nside_out ${band_name} $beam_arcmin_sidelobe -${i} n $unit y $dir/$filename_synch $dir/$filename_dust
done


# CMB+Synch+Dust+Dipole
filename_cmb=cmb_map_$band
filename_dust=thermaldust_map_$band
filename_synch=synchrotron_map_$band
for((i=0 ; i<=100 ; i+=10))
do 
    filename_out=${filename_cmb}_FG_Dipole_sidelobe_${i}dB_b${beam_arcmin_sidelobe}arcmin_ud${nside_out}
    bsub -o $dir_h/log/log19_${i}_${band_name}_${unit}.txt -q e python gen_sidelobemap_dipole.py $dir/$filename_cmb $dir_out/$filename_out $nside_out ${band_name} $beam_arcmin_sidelobe -${i} y $unit y $dir/$filename_synch $dir/$filename_dust
done

