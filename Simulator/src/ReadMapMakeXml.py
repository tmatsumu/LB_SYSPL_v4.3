'''
this was originally ReadScanXml_ver2.py
'''
from xml.dom import minidom, Node
from numpy import *
import sys

def convert_Gregorian2Julian(year, month, day, hour, minute, sec):
    a = (14-month)/12
    a = int(a)
    y = year + 4800 - a
    m = month + 12*a - 3
    JDN = day+int((153*m+2)/5)+365*y+int(y/4)-int(y/100)+int(y/400)-32045
    JD = JDN + (hour-12)/24. + minute/1440. + sec/86400.
    return array(JD)

def readxml_date(node):
    for child in node.childNodes:
        if child.nodeType == Node.ELEMENT_NODE:
            if child.tagName == 'date_i': date_i = getText(child)
            if child.tagName == 'date_f': date_f = getText(child)
    y_i = date_i[0:4];    m_i = date_i[4:6];    d_i = date_i[6:]
    y_f = date_f[0:4];    m_f = date_f[4:6];    d_f = date_f[6:]
    MJD_i = convert_Gregorian2Julian(int(y_i), int(m_i), int(d_i), 0.,0.,0.)-2400000.5
    MJD_f = convert_Gregorian2Julian(int(y_f), int(m_f), int(d_f), 0.,0.,0.)-2400000.5
    return MJD_i, MJD_f

# <log>
def readxml_log(node):
    for child in node.childNodes:
        if child.nodeType == Node.ELEMENT_NODE:
            if child.tagName == 'runID': runID = getText(child)
            if child.tagName == 'dir_simulator': dir_simulator = getText(child)
            if child.tagName == 'debug': debug = getText(child)
            if child.tagName == 'machine': machine = getText(child)
            if child.tagName == 'silent': silent = getText(child)
    return runID, dir_simulator, debug, machine, silent
# </log>

# <basicpar>
def readxml_basicpar(node):
    for child in node.childNodes:
        if child.nodeType == Node.ELEMENT_NODE:
            if child.tagName == 'nside': nside = getText(child)
            if child.tagName == 'run_type': run_type = getText(child)
            if child.tagName == 'coord': coord = getText(child)
            if child.tagName == 'pixelmapio': pixelmapio = getText(child)
            if child.tagName == 'gen_tod': gen_tod = getText(child)
            if child.tagName == 'TQU': TQU = getText(child)
    return int(nside), run_type, coord, pixelmapio, gen_tod, TQU
# </basicpar>

# <simulations>
def readxml_simulations(node,TQU):
    for child in node.childNodes:
        if child.nodeType == Node.ELEMENT_NODE:
            if ((TQU == 'T') or (TQU=='TQU')):
                if child.tagName == 'file_input_maps': file_input_maps = getText(child)
                if child.tagName == 'file_input_noise': file_input_noise = getText(child)
                if child.tagName == 'file_input_fpdb': file_input_fpdb = getText(child)
                if child.tagName == 'file_input_muellermatrix': file_input_muellermatrix = getText(child)
            if TQU == 'TQU2':
                if child.tagName == 'file_input_maps': file_input_maps = getText(child)
                if child.tagName == 'file_input_maps2': file_input_maps2 = getText(child)
                if child.tagName == 'file_input_noise': file_input_noise = getText(child)
                if child.tagName == 'file_input_fpdb': file_input_fpdb = getText(child)
                if child.tagName == 'file_input_muellermatrix': file_input_muellermatrix = getText(child)
            if TQU == 'T_bandpassmismatch':
                if child.tagName == 'file_input_maps': file_input_maps = getText(child)
                if child.tagName == 'file_input_map_dust': file_input_map_dust = getText(child)
                if child.tagName == 'file_input_map_synch': file_input_map_synch = getText(child)
                if child.tagName == 'file_input_bandpassmismatch': file_input_bandpassmismatch = getText(child)
                if child.tagName == 'file_input_noise': file_input_noise = getText(child); print file_input_noise
                if child.tagName == 'file_input_fpdb': file_input_fpdb = getText(child)
                if child.tagName == 'file_input_muellermatrix': file_input_muellermatrix = getText(child)

    if ((TQU == 'T') or (TQU=='TQU')):
        return file_input_maps, file_input_noise, file_input_fpdb, file_input_muellermatrix
    if TQU == 'TQU2':    
        return file_input_maps, file_input_maps2, file_input_noise, file_input_fpdb, file_input_muellermatrix
    if TQU == 'T_bandpassmismatch':
        return file_input_maps, file_input_noise, file_input_fpdb, file_input_muellermatrix, file_input_map_dust, file_input_map_synch, file_input_bandpassmismatch

# </simulations>

# <database>
def readxml_database(node):
    for child in node.childNodes:
        if child.nodeType == Node.ELEMENT_NODE:
            if child.tagName == 'file_fpdb_mmin': file_fpdb_mmin = getText(child)
            if child.tagName == 'file_muellermatrix': file_muellermatrix = getText(child)
#            if child.tagName == 'file_relgain': file_relgain = getText(child)
#            if child.tagName == 'file_boloid': file_boloid = getText(child)
#            if child.tagName == 'file_flag_pixel': file_flag_pixel = getText(child)
            if child.tagName == 'db_ces': db_ces = getText(child)
            if child.tagName == 'db_ces2': db_ces2 = getText(child)
            if child.tagName == 'db_gain': db_gain = getText(child)
            if child.tagName == 'gain_type': gain_type = getText(child)
            if child.tagName == 'gain_corr': gain_corr = getText(child)
#    return file_fpdb_mmin, file_muellermatrix, file_relgain, file_boloid, file_flag_pixel, db_ces, db_ces2, db_gain, gain_type
    return file_fpdb_mmin, file_muellermatrix, db_ces, db_ces2, db_gain, gain_type, gain_corr
# </database>

# <filetering_choice>
def readxml_filtering_choice(node):
    for child in node.childNodes:
        if child.nodeType == Node.ELEMENT_NODE:
            if child.tagName == 'filter_choice': filter_choice = getText(child)
            if child.tagName == 'poly': poly = getText(child)
#            if child.tagName == 'file_noisefft': file_noisefft = getText(child)
    return filter_choice, int(poly)
# </filetering_choice>

# <output_file>
def readxml_output_file(node):
    for child in node.childNodes:
        if child.nodeType == Node.ELEMENT_NODE:
            if child.tagName == 'dir_simedmap': dir_simedmap = getText(child)
            if child.tagName == 'dir_combinedmap': dir_combinedmap = getText(child)
            if child.tagName == 'Tn_map': Tn_map = getText(child)
            if child.tagName == 'Td_map': Td_map = getText(child)
            if child.tagName == 'AA_map': AA_map = getText(child)
            if child.tagName == 'AB_map': AB_map = getText(child)
            if child.tagName == 'BB_map': BB_map = getText(child)
            if child.tagName == 'Ad_map': Ad_map = getText(child)
            if child.tagName == 'Bd_map': Bd_map = getText(child)
    return dir_simedmap, dir_combinedmap, Tn_map, Td_map, AA_map, AB_map, BB_map, Ad_map, Bd_map
# </output_file>

def getText(node):
    s = ''
    for child in node.childNodes:
        if child.nodeType == Node.TEXT_NODE:
            s += child.wholeText
    return s

def Get_Mapmake_Inputs(filename):
    print '[Get_Mapmake_Inputs in ReadMapXml.py]', filename
    doc = minidom.parse(filename)

    for node in doc.getElementsByTagName('log'):
        runID, dir_simulator, debug, machine, silent = readxml_log(node)
    for node in doc.getElementsByTagName('basicpar'):
        nside, run_type, coord, pixelmapio, gen_tod, TQU = readxml_basicpar(node)
    for node in doc.getElementsByTagName('simulations'):
        if ((TQU=='T') or (TQU=='TQU')):
            file_input_maps, file_input_noise, file_input_fpdb, file_input_muellermatrix = readxml_simulations(node, TQU)
        if TQU=='TQU2':
            file_input_maps, file_input_maps2, file_input_noise, file_input_fpdb, file_input_muellermatrix = readxml_simulations(node, TQU)
        if TQU=='T_bandpassmismatch':
            file_input_maps, file_input_noise, file_input_fpdb, file_input_muellermatrix, file_input_dust, file_input_synch, file_input_bandpassmismatch  = readxml_simulations(node, TQU)
    for node in doc.getElementsByTagName('database'):
        file_fpdb_mmin, file_muellermatrix, db_ces, db_ces2, db_gain, gain_type, gain_corr = readxml_database(node)
#        file_fpdb_mmin, file_muellermatrix, file_relgain, file_boloid, file_flag_pixel, db_ces, db_ces2, db_gain, gain_type = readxml_database(node)
    for node in doc.getElementsByTagName('filtering_choice'):
        filter_choice, poly = readxml_filtering_choice(node)
#        choice, poly, file_noisefft = readxml_filtering_choice(node)
    for node in doc.getElementsByTagName('output_file'):
        dir_simedmap, dir_combinedmap, Tn, Td, AA, BB, AB, Ad, Bd = readxml_output_file(node)

    if ((TQU=='T') or (TQU=='TQU')):
        out = {'runID': runID, 
               'dir_simulator': dir_simulator,
               'debug': debug,
               'machine': machine,
               'silent': silent,
               'nside':int(nside), 'run_type':run_type, 'coord':coord, 
               'pixelmapio': pixelmapio, 'gen_tod': gen_tod,
               'file_input_maps':file_input_maps, 
               'file_input_noise':file_input_noise,
               'file_input_fpdb':file_input_fpdb,
               'file_input_muellermatrix':file_input_muellermatrix,
               'TQU':TQU,
               'file_fpdb_mmin':file_fpdb_mmin,
               'file_muellermatrix':file_muellermatrix,
               'filter_choice': filter_choice,
               'poly':poly,
               'db_ces':db_ces,
               'db_ces2':db_ces2,
               'db_gain':db_gain,
               'gain_type':gain_type,
               'gain_corr':gain_corr,
               'dir_simedmap':dir_simedmap, 
               'dir_combinedmap':dir_combinedmap, 
               'fname_Tn':Td, 'fname_Td':Td, 
               'fname_AA':AA, 'fname_BB':BB, 'fname_AB':AB, 'fname_Ad':Ad, 'fname_Bd':Bd}
        

    if TQU=='TQU2':
        out = {'runID': runID, 
               'dir_simulator': dir_simulator,
               'debug': debug,
               'machine': machine,
               'silent': silent,
               'nside':int(nside), 'run_type':run_type, 'coord':coord, 
               'pixelmapio': pixelmapio, 'gen_tod': gen_tod,
               'file_input_maps':file_input_maps, 
               'file_input_maps2':file_input_maps2, 
               'file_input_noise':file_input_noise,
               'file_input_fpdb':file_input_fpdb,
               'file_input_muellermatrix':file_input_muellermatrix,
               'TQU':TQU,
               'file_fpdb_mmin':file_fpdb_mmin,
               'file_muellermatrix':file_muellermatrix,
               'filter_choice': filter_choice,
               'poly':poly,
               'db_ces':db_ces,
               'db_ces2':db_ces2,
               'db_gain':db_gain,
               'gain_type':gain_type,
               'gain_corr':gain_corr,
               'dir_simedmap':dir_simedmap, 
               'dir_combinedmap':dir_combinedmap, 
               'fname_Tn':Td, 'fname_Td':Td, 
               'fname_AA':AA, 'fname_BB':BB, 'fname_AB':AB, 'fname_Ad':Ad, 'fname_Bd':Bd}
        

    if TQU=='T_bandpassmismatch':
        out = {'runID': runID, 
               'dir_simulator': dir_simulator,
               'debug': debug,
               'machine': machine,
               'silent': silent,
               'nside':int(nside), 'run_type':run_type, 'coord':coord, 
               'pixelmapio': pixelmapio, 'gen_tod': gen_tod,
               'file_input_maps':file_input_maps, 
               'file_input_dust':file_input_dust, 
               'file_input_synch':file_input_synch, 
               'file_input_noise':file_input_noise,
               'file_input_bandpassmismatch':file_input_bandpassmismatch,
               'file_input_fpdb':file_input_fpdb,
               'file_input_muellermatrix':file_input_muellermatrix,
               'TQU':TQU,
               'file_fpdb_mmin':file_fpdb_mmin,
               'file_muellermatrix':file_muellermatrix,
               'filter_choice': filter_choice,
               'poly':poly,
               'db_ces':db_ces,
               'db_ces2':db_ces2,
               'db_gain':db_gain,
               'gain_type':gain_type,
               'gain_corr':gain_corr,
               'dir_simedmap':dir_simedmap, 
               'dir_combinedmap':dir_combinedmap, 
               'fname_Tn':Td, 'fname_Td':Td, 
               'fname_AA':AA, 'fname_BB':BB, 'fname_AB':AB, 'fname_Ad':Ad, 'fname_Bd':Bd}
        
    return out

if __name__ == '__main__': _main()


