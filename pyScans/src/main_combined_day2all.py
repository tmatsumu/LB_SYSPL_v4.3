import numpy as np
import pylab as py
import sqlite3 as sq
import os.path
import sys

pi = np.pi
radeg = (180./pi)

def read_lb_observation(filename,sq_command):
    print filename
    conn = sq.connect(filename)
    c = conn.cursor()
#    c.execute('.table;')
#    print c
#    c.close()
#    sys.exit()
    c.execute(sq_command)
    id=[];juliantime=[];dir_ptg=[];
    for ar in c:
        id.append(int(ar[0]))
        juliantime.append(float(ar[1]))
        dir_ptg.append(ar[2])
    c.close()
    CESdb = {'id':id,'juliantime':juliantime,'dir_ptg':dir_ptg}
    return CESdb

def display_all(db_dict):
    keys = db_dict.keys()
    num = len(db_dict[keys[0]])
    print keys
    for i in range(num):
        tmp = []
        for j in keys:
            tmp.append(db_dict[j][i])
            print db_dict[j][i]

# read the pointing file                                                                                                                               
#  input: boresight fits filename                                                                                                                      
#  output: {ra, dec, pa, hwp}                                                                                                                          
def read_LBptg_writeinRADEC(filename):
    print ""
    print "[READ PTG]: BEGIN reading LB pointing: "+filename
    ptg = np.load(filename)
    lat = ptg['lat']
    lon = ptg['lon']
    pa = ptg['pa']
    hwp = np.zeros(len(lat))
    ra,dec = euler_astrolib(lon*radeg,lat*radeg,4,FK4='J2000')
    ptg_package = {'pa':pa, 'ra':ra/radeg, 'dec':dec/radeg, 'hwp':hwp}
#    print 'HWP angle [degs]', ptg[4]/pi*180.                                                                                                          
    print "[READ PTG]: END reading "
    return ptg_package


def euler_astrolib(ai,bi,select, FK4=False, radian=False):
    '''                                                                                                                                          
    ; NAME:                                                                                                                                      
    ;     EULER                                                                                                                                  
    ; PURPOSE:                                                                                                                                   
    ;     Transform between Galactic, celestial, and ecliptic coordinates.                                                                       
    ; EXPLANATION:                                                                                                                               
    ;     Use the procedure ASTRO to use this routine interactively                                                                              
    ;                                                                                                                                            
    ; CALLING SEQUENCE:                                                                                                                          
    ;      EULER, AI, BI, AO, BO, [ SELECT, /FK4 ]                                                                                               
    ;                                                                                                                                            
    ; INPUTS:                                                                                                                                    
    ;       AI - Input Longitude in DEGREES, scalar or vector.  If only two                                                                      
    ;               parameters are supplied, then  AI and BI will be modified to                                                                 
    ;               contain the output longitude and latitude.                                                                                   
    ;       BI - Input Latitude in DEGREES                                                                                                       
    ;                                                                                                                                            
    ; OPTIONAL INPUT:                                                                                                                            
    ;       SELECT - Integer (1-6) specifying type of coordinate transformation.                                                                 
    ;                                                                                                                                            
    ;      SELECT   From          To        |   SELECT      From            To                                                                   
    ;       1     RA-Dec (2000)  Galactic   |     4       Ecliptic      RA-Dec                                                                   
    ;       2     Galactic       RA-DEC     |     5       Ecliptic      Galactic                                                                 
    ;       3     RA-Dec         Ecliptic   |     6       Galactic      Ecliptic                                                                 
    ;                                                                                                                                            
    ;      If omitted, program will prompt for the value of SELECT                                                                               
    ;      Celestial coordinates (RA, Dec) should be given in equinox J2000                                                                      
    ;      unless the /FK4 keyword is set.                                                                                                       
    ; OUTPUTS:                                                                                                                                   
    ;       AO - Output Longitude in DEGREES                                                                                                     
    ;       BO - Output Latitude in DEGREES                                                                                                      
    ;                                                                                                                                            
    ; INPUT KEYWORD:                                                                                                                             
    ;       /FK4 - If this keyword is set and non-zero, then input and output                                                                    
    ;             celestial and ecliptic coordinates should be given in equinox                                                                  
    ;             B1950.                                                                                                                         
    ;                                                                                                                                            
    ; NOTES:                                                                                                                                     
    ;       EULER was changed in December 1998 to use J2000 coordinates as the                                                                   
    ;       default, ** and may be incompatible with earlier versions***.                                                                        
    ; REVISION HISTORY:                                                                                                                          
    ;       Written W. Landsman,  February 1987                                                                                                  
    ;       Adapted from Fortran by Daryl Yentis NRL                                                                                             
    ;       Converted to IDL V5.0   W. Landsman   September 1997                                                                                 
    ;       Made J2000 the default, added /FK4 keyword  W. Landsman December 1998                                                                
    ;-                                                                                                                                           
    '''
    twopi   =   2.0*pi
    fourpi  =   4.0*pi
    deg_to_rad = (180.0/pi)

#;   J2000 coordinate conversions are based on the following constants                                                                           
#;  eps = 23.4392911111d              Obliquity of the ecliptic                                                                                  
#;  alphaG = 192.85948d               Right Ascension of Galactic North Pole                                                                     
#;  deltaG = 27.12825d                Declination of Galactic North Pole                                                                         
#;  lomega = 32.93192d                Galactic longitude of celestial equator                                                                    
#;  alphaE = 180.02322d              Ecliptic longitude of Galactic North Pole                                                                   
#;  deltaE = 29.811438523d            Ecliptic latitude of Galactic North Pole                                                                   
#;  Eomega  = 6.3839743d              Galactic longitude of ecliptic equator                                                                     

    if FK4:
        equinox = '(B1950)'
        psi   = np.array([ 0.57595865315, 4.9261918136, 0.00000000000, 0.0000000000, 0.11129056012, 4.7005372834])
        stheta = np.array([ 0.88781538514,-0.88781538514, 0.39788119938,-0.39788119938, 0.86766174755,-0.86766174755])
        ctheta = np.array([ 0.46019978478, 0.46019978478, 0.91743694670, 0.91743694670, 0.49715499774, 0.49715499774])
        phi  = np.array([ 4.9261918136,  0.57595865315, 0.0000000000, 0.00000000000, 4.7005372834, 0.11129056012])

    equinox = '(J2000)'
    psi   =  np.array([ 0.57477043300, 4.9368292465, 0.00000000000, 0.0000000000, 0.11142137093, 4.71279419371])
    stheta = np.array([ 0.88998808748,-0.88998808748, 0.39777715593,-0.39777715593, 0.86766622025,-0.86766622025])
    ctheta = np.array([ 0.45598377618, 0.45598377618, 0.91748206207, 0.91748206207, 0.49714719172, 0.49714719172])
    phi  = np.array([ 4.9368292465,  0.57477043300, 0.0000000000, 0.00000000000, 4.71279419371, 0.11142137093])

    i  = select - 1                         # IDL offset                                                                                         
    a  = ai/deg_to_rad - phi[i]
    b = bi/deg_to_rad
    sb = np.sin(b)
    cb = np.cos(b)
    cbsa = cb * np.sin(a)
    b  = -stheta[i] * cbsa + ctheta[i] * sb
    bo = np.arcsin(b)*deg_to_rad

    a =  np.arctan2( ctheta[i] * cbsa + stheta[i] * sb, cb * np.cos(a) )
    ao = ( (a+psi[i]+fourpi) % twopi) * deg_to_rad

    return ao, bo

def write_pointing_biary(filename_out,data):
    num = len(data)
    f = open(filename_out, "ab")
    for i in range(0, num) :
        f.write(data[i])
    f.close

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#dir_db  = sys.argv[1]
#sub_dir = sys.argv[2]
filename_db = sys.argv[1]
sq_command = sys.argv[2]

filename_out = sys.argv[3]

DB_arr= read_lb_observation(filename_db+'.db', sq_command)

display_all(DB_arr)

if os.path.isfile(filename_out+'_ra.bi'):
    print 'The output file already exists!'
    print filename_out+'_ra.bi'
    sys.exit()
if os.path.isfile(filename_out+'_dec.bi'):
    print 'The output file already exists!'
    print filename_out+'_dec.bi'
    sys.exit()
if os.path.isfile(filename_out+'_psi.bi'):
    print 'The output file already exists!'
    print filename_out+'_psi.bi'
    sys.exit()

nb_scanset = len(DB_arr['id'])
for i_scanset in range(0,nb_scanset):
    print DB_arr['dir_ptg'][i_scanset]
    ptg_package = read_LBptg_writeinRADEC(DB_arr['dir_ptg'][i_scanset]+'.npz')
    write_pointing_biary(filename_out+'_ra.bi', ptg_package['ra'])
    write_pointing_biary(filename_out+'_dec.bi', ptg_package['dec'])
    write_pointing_biary(filename_out+'_psi.bi', ptg_package['pa'])
