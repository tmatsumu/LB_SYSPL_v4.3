#!/bin/sh

runID="test"

#----------------------------------------------------
dir_bin=/raid/users/tmatsumu/work/develop/PBI/PB1_NTP/src
dir_xml=/raid/users/tmatsumu/work/develop/PBI/PB1_NTP/xml

#----------------------------------------------------

cat $dir_xml/header_template.xml   > $dir_xml/xml_${runID}.xml
echo '<Mapmaking>'                  >> $dir_xml/xml_${runID}.xml
cat $dir_xml/Global_set.xml       >> $dir_xml/xml_${runID}.xml

cat $dir_xml/misc.xml             >> $dir_xml/xml_${runID}.xml
#echo '<data_selection>'           >> $dir_xml/xml_${runID}.xml
#python $dir_bin/selected_data.py  >> $dir_xml/xml_${runID}.xml
#echo '</data_selection>'          >> $dir_xml/xml_${runID}.xml
#
#echo '<bolo_selection>'           >> $dir_xml/xml_${runID}.xml
#python $dir_bin/selected_bolo.py  >> $dir_xml/xml_${runID}.xml 
#echo '</bolo_selection>'          >> $dir_xml/xml_${runID}.xml
echo ''                             >> $dir_xml/xml_${runID}.xml
echo '</Mapmaking>'                 >> $dir_xml/xml_${runID}.xml
#----------------------------------------------------

python $dir_bin/main_mapmaker_dist.py $dir_xml/xml_${runID}.xml test