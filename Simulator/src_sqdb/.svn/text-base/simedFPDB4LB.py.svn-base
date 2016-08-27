import numpy as np
import pylab as py
import fileinput
import sys

'''
simedFPDB4LB.py 
this code was originally copied from simedFPDB_ver3.py that was written for PB.
T. Matsumura, 2013-3-30
'''

pi = np.pi

class move():
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.num = len(self.x)

    def trans_rot(self, angle, dist_x, dist_y):
        self.angle = angle 
        xx, yy = np.zeros(self.num), np.zeros(self.num)
        for i in range(self.num):
            xx[i] = self.x[i]*np.cos(self.angle)-self.y[i]*np.sin(self.angle)
            yy[i] = self.x[i]*np.sin(self.angle)+self.y[i]*np.cos(self.angle)
        return xx+dist_x, yy+dist_y

    def flip(self):
        self.az_f  = - self.x
        self.el_f  = - self.y
        return self.az_f, self.el_f
    
def label():
#    tag = ['10.2','10.4','8.2.0','10.1','9.4','10.5','10.3']
    tag = ['HFW1','HFW2','HFW3','HFW4','HFW5']

    wedge = []
    idx = []
    tag_out = []
    pol_state = []
    for i in range(0,7):
        for j in range(0,91):
            tag_tmp_t = tag[i]+'_'+str(j)+'t'
            tag_tmp_b = tag[i]+'_'+str(j)+'b'
            tag_out.append(tag_tmp_t)
            tag_out.append(tag_tmp_b)
            idx.append(i*j+j+1)
            wedge.append(tag[i])
            pol_state_t = 't'
            pol_state_b = 'b'
            pol_state.append(pol_state_t)
            pol_state.append(pol_state_b)
    return tag_out, idx, wedge, pol_state

def read_FPDB_designed(filename):
    wedge= []
    bolo = []
    pix = []
    torb = []
    polang = []
    x_w = []
    y_w = []
    x_b = []
    y_b = []
    r_b = []
    theta_x = []
    theta_z = []
        
    i = 0
    for line in fileinput.input(filename):
        ar = line.split()
        if ((len(ar)>1) & (i>0)):
            wedge.append(int(ar[0]))
            bolo.append(int(ar[1]))
            pix.append(int(ar[2]))
            torb.append(ar[3])
            polang.append(float(ar[4]))
            x_w.append(float(ar[5]))
            y_w.append(float(ar[6]))
            x_b.append(float(ar[7]))
            y_b.append(float(ar[8]))
            r_b.append(float(ar[9]))
            theta_x.append(float(ar[10]))
            theta_z.append(float(ar[11]))
        i += 1
    print "[READ FPDB]: End reading "+filename
    FPDBlist = {'wedge#': np.array(wedge), 'bolo': np.array(bolo), 'pix': np.array(pix), 'torb': np.array(torb),
                'polang': np.array(polang), 'x_w': np.array(x_w), 'y_w': np.array(y_w),
                'x_b': np.array(x_b), 'y_b': np.array(y_b),  'r_b': np.array(r_b),
                'theta_x': np.array(theta_x), 'theta_z': np.array(theta_z)}
    return FPDBlist
 
def gen_polvec(az,el,polang,amp):
    az_out1 = az+amp*np.cos(polang)
    az_out2 = az-amp*np.cos(polang)
    el_out1 = el+amp*np.sin(polang)
    el_out2 = el-amp*np.sin(polang)
    return np.array([az_out1,az_out2]), np.array([el_out1,el_out2])

def write_txt(filenameout,x0,y0,polang0,x1,y1,polang1,x2,y2,polang2,x3,y3,polang3,x4,y4,polang4):
    x = np.hstack((x0,x1,x2,x3,x4))
    y = np.hstack((y0,y1,y2,y3,y4))
    polang = np.hstack((polang0,polang1,polang2,polang3,polang4))
    num = len(x)
#    num_wpix = 91
#    num_wbolo = num_wpix *2 
#    num_w = 7
    num_wpix = 37
    num_wbolo = num_wpix *2 
    num_w = 5
    num_tot = num_wpix*num_w

    wafer = np.array(['HFW1','HFW2','HFW3','HFW4','HFW5'])
#    wafer = np.array(['10.2','10.4','8.2.0','10.1','9.4','10.5','10.3'])

    j = 0
    k = 0

    f = open(filenameout, "w")    
    f.write('index pixel# az [degs] el [degs] polang [degs] wafer boloID \n')
    for i in range(0,num_tot):
        f.write('%d   %d   %f   %f   %4.1f   %s   %s \n' % (2*i, k+1, x[i], y[i], polang[2*i], wafer[j], wafer[j]+'_'+str(k+1)+'t'))
        f.write('%d   %d   %f   %f   %4.1f   %s   %s \n' % (2*i+1, k+1, x[i], y[i], polang[2*i+1], wafer[j], wafer[j]+'_'+str(k+1)+'b'))
        k=k+1
        if (i+1)%91 == 0: 
            j=j+1
            k=0
    f.close()

def write_txt_tmp(filenameout,x0,y0,polang0,x1,y1,polang1,x2,y2,polang2,x3,y3,polang3,x4,y4,polang4):
    x = np.hstack((x0,x1,x2,x3,x4))
    y = np.hstack((y0,y1,y2,y3,y4))
    polang = np.hstack((polang0,polang1,polang2,polang3,polang4))
    num = len(x)
#    num_wpix = 91
#    num_wbolo = num_wpix *2 
#    num_w = 7
    num_wpix = 37
    num_wbolo = num_wpix *2 
    num_w = 5
    num_tot = num_wpix*num_w

#    wafer = np.array(['10.2','10.4','8.2.0','10.1','9.4','10.5','10.3'])
    wafer = np.array(['HFW1','HFW2','HFW3','HFW4','HFW5'])

    j = 0
    k = 0

    beam = 3.5/60.
    amp = 1.
    poleff = 1.
    f = open(filenameout, "w")    
#    f.write('index pixel# az [degs] el [degs] polang [degs] wafer boloID \n')
    f.write('ch, az[deg], el[deg], sig_x[deg], sig_y[deg], amp, polang[deg], poleff, wafer, pix, torb, flag \n')
    for i in range(0,num_tot):
        f.write('%d   %f   %f   %f   %f   %1.1f   %4.1f  %f  %s %d  %s %d\n' % (2*i, x[i], y[i], beam, beam, amp, polang[2*i], poleff, wafer[j], k+1, 't', 1))
        f.write('%d   %f   %f   %f   %f   %1.1f   %4.1f  %f  %s %d  %s %d\n' % (2*i+1, x[i], y[i], beam, beam, amp, polang[2*i+1], poleff, wafer[j], k+1, 'b', 1))
        k=k+1
        if (i+1)%91 == 0: 
            j=j+1
            k=0
    f.close()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def run():
    skyview=True
#    skyview=False
    
    if skyview==True:
        plate_scale = (1./1100.) * (180./pi) # [rad/mm] -> [deg/mm]

    if skyview==False:
        plate_scale = 1.


    pix_sep = 12.*plate_scale
    az_sep = pix_sep
    el_sep = pix_sep*np.sqrt(3.)/2.
    
    sep_dist= 85.*plate_scale
    if skyview==True: pol_vec_length = 2.*plate_scale
    if skyview==False: pol_vec_length = 2.

#+++++++++++++++++++++

#    columns = np.array([6,7,8,9,10,11,10,9,8,7,6])
    columns = np.array([4,5,6,7,6,5,4])
    columns_i = -az_sep*(np.float_(columns)-1.)/2.
    num_columns = len(columns)
    
    rows = np.array([3,2,1,0,-1,-2,-3])
    num_rows = len(rows)
    
    # type1: pix1 (t=-90,b=0)
    # type-1: pix2 (t=-135,b=-45)
    polang_pixindex = np.array([1,5,10,16,23,29,34])
    polang_index = np.array([1,-1,1,-1,-1,-1,-1])
    
    az = []
    el = []
    idx = []
    polang = []
    polang_t = []
    polang_b = []
    
    polang_index_tmp = 1
    idx_tmp = 1
    for j in range(0,num_rows):
        for i in range(0,columns[j]):
            #        print j, i, idx_tmp, columns_i[j], columns_i[j]+az_sep*i, rows[j], polang_index_tmp
            idx.append(idx_tmp)
            az_i = (columns_i[j]+az_sep*i)
            az.append(az_i)
            el.append(el_sep*rows[j])

            if ((idx_tmp == polang_pixindex[j]) and (polang_index[j] == 1)): 
                pol_t = -90./180.*pi
                pol_b = 0./180.*pi
                polang.append(pol_t)
                polang.append(pol_b)
                polang_t.append(pol_t)
                polang_b.append(pol_b)
                polang_index_tmp = 1

            if ((idx_tmp != polang_pixindex[j]) and (polang_index_tmp == 1)):
                pol_t = -90./180.*pi
                pol_b = 0./180.*pi
                polang.append(pol_t)
                polang.append(pol_b)
                polang_t.append(pol_t)
                polang_b.append(pol_b)
                polang_index_tmp = 1

            if ((idx_tmp == polang_pixindex[j]) and (polang_index[j] == -1)):
                pol_t = -135./180.*pi
                pol_b = -45./180.*pi
                polang.append(pol_t)
                polang.append(pol_b)
                polang_t.append(pol_t)
                polang_b.append(pol_b)
                polang_index_tmp = -1

            if ((idx_tmp != polang_pixindex[j]) and (polang_index_tmp == -1)):
                pol_t = -135./180.*pi
                pol_b = -45./180.*pi
                polang.append(pol_t)
                polang.append(pol_b)
                polang_t.append(pol_t)
                polang_b.append(pol_b)
                polang_index_tmp = -1

            polang_index_tmp = polang_index_tmp*(-1)
            idx_tmp += 1

    az = az[::-1]

# center wafer, wafer#-0
    out0 = move(az,el)
    wafer_angle = 90./180.*pi
    x0, y0 = out0.trans_rot(wafer_angle,0.,0.)
    polvec_x0_t, polvec_y0_t = gen_polvec(x0,y0,np.float_(polang_t)+wafer_angle,pol_vec_length)
    polvec_x0_b, polvec_y0_b = gen_polvec(x0,y0,np.float_(polang_b)+wafer_angle,pol_vec_length)
    polang0 = np.float_(polang)+wafer_angle/pi*180.

# top wafer in the photon view, wafer#-1 
#    wafer_angle = 30./180.*pi
#    wafer_pos = 120./180.*pi
    wafer_angle = -30./180.*pi
    wafer_pos = 60./180.*pi
    out1 = move(az,el)
    x1, y1 = out1.trans_rot(wafer_angle,sep_dist*np.cos(wafer_pos),sep_dist*np.sin(wafer_pos))
    polvec_x1_t, polvec_y1_t = gen_polvec(x1,y1,np.float_(polang_t)+wafer_angle,pol_vec_length)
    polvec_x1_b, polvec_y1_b = gen_polvec(x1,y1,np.float_(polang_b)+wafer_angle,pol_vec_length)
    polang1 = np.float_(polang)+wafer_angle/pi*180.

    wafer_angle = -90./180.*pi
    wafer_pos = 0./180.*pi
    out2 = move(az,el)
    x2, y2 = out2.trans_rot(wafer_angle,sep_dist*np.cos(wafer_pos),sep_dist*np.sin(wafer_pos))
    polvec_x2_t, polvec_y2_t = gen_polvec(x2,y2,np.float_(polang_t)+wafer_angle,pol_vec_length)
    polvec_x2_b, polvec_y2_b = gen_polvec(x2,y2,np.float_(polang_b)+wafer_angle,pol_vec_length)
    polang2 = np.float_(polang)+wafer_angle/pi*180.

    wafer_angle = 210./180.*pi
    wafer_pos = -60./180.*pi
    out3 = move(az,el)
    x3, y3 = out3.trans_rot(wafer_angle,sep_dist*np.cos(wafer_pos),sep_dist*np.sin(wafer_pos))
    polvec_x3_t, polvec_y3_t = gen_polvec(x3,y3,np.float_(polang_t)+wafer_angle,pol_vec_length)
    polvec_x3_b, polvec_y3_b = gen_polvec(x3,y3,np.float_(polang_b)+wafer_angle,pol_vec_length)
    polang3 = np.float_(polang)+wafer_angle/pi*180.

    wafer_angle = 150./180.*pi
    wafer_pos = -120./180.*pi
    out4 = move(az,el)
    x4, y4 = out4.trans_rot(wafer_angle,sep_dist*np.cos(wafer_pos),sep_dist*np.sin(wafer_pos))
    polvec_x4_t, polvec_y4_t = gen_polvec(x4,y4,np.float_(polang_t)+wafer_angle,pol_vec_length)
    polvec_x4_b, polvec_y4_b = gen_polvec(x4,y4,np.float_(polang_b)+wafer_angle,pol_vec_length)
    polang4 = np.float_(polang)+wafer_angle/pi*180.

#    wafer_angle = 120./180.*pi
#    wafer_pos = -150./180.*pi
#    out5 = move(az,el)
#    x5, y5 = out5.trans_rot(wafer_angle,sep_dist*np.cos(wafer_pos),sep_dist*np.sin(wafer_pos))
#    polvec_x5_t, polvec_y5_t = gen_polvec(x5,y5,np.float_(polang_t),pol_vec_length)
#    polvec_x5_b, polvec_y5_b = gen_polvec(x5,y5,np.float_(polang_b),pol_vec_length)
#    polang5 = np.float_(polang)+wafer_angle/pi*180.

#    wafer_angle = 60./180.*pi
#    wafer_pos = 150./180.*pi
#    out6 = move(az,el)
#    x6, y6 = out6.trans_rot(wafer_angle,sep_dist*np.cos(wafer_pos),sep_dist*np.sin(wafer_pos))
#    polvec_x6_t, polvec_y6_t = gen_polvec(x6,y6,np.float_(polang_t),pol_vec_length)
#    polvec_x6_b, polvec_y6_b = gen_polvec(x6,y6,np.float_(polang_b),pol_vec_length)
#    polang6 = np.float_(polang)+wafer_angle/pi*180.

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if skyview==True:
        out = move(x0,y0); x0,y0=out.flip()
        out = move(x1,y1); x1,y1=out.flip()
        out = move(x2,y2); x2,y2=out.flip()
        out = move(x3,y3); x3,y3=out.flip()
        out = move(x4,y4); x4,y4=out.flip()

        out = move(polvec_x0_t,polvec_y0_t); polvec_x0_t,polvec_y0_t=out.flip()
        out = move(polvec_x0_b,polvec_y0_b); polvec_x0_b,polvec_y0_b=out.flip()
        out = move(polvec_x1_t,polvec_y1_t); polvec_x1_t,polvec_y1_t=out.flip()
        out = move(polvec_x1_b,polvec_y1_b); polvec_x1_b,polvec_y1_b=out.flip()
        out = move(polvec_x2_t,polvec_y2_t); polvec_x2_t,polvec_y2_t=out.flip()
        out = move(polvec_x2_b,polvec_y2_b); polvec_x2_b,polvec_y2_b=out.flip()
        out = move(polvec_x3_t,polvec_y3_t); polvec_x3_t,polvec_y3_t=out.flip()
        out = move(polvec_x3_b,polvec_y3_b); polvec_x3_b,polvec_y3_b=out.flip()
        out = move(polvec_x4_t,polvec_y4_t); polvec_x4_t,polvec_y4_t=out.flip()
        out = move(polvec_x4_b,polvec_y4_b); polvec_x4_b,polvec_y4_b=out.flip()

#        out = move(x5,y5); x5,y5=out.flip()
#        out = move(x6,y6); x6,y6=out.flip()
#        write_txt('simfullfpdb_skview.txt',x0,y0,polang0,x1,y1,polang1,x2,y2,polang2,x3,y3,polang3,x4,y4,polang4,x5,y5,polang5,x6,y6,polang6)
#        write_txt_tmp('FPDB4sim_simfullfpdb_skyview.txt',x0,y0,polang0,x1,y1,polang1,x2,y2,polang2,x3,y3,polang3,x4,y4,polang4,x5,y5,polang5,x6,y6,polang6)

        write_txt('simfullfpdb_skview.txt',x0,y0,polang0,x1,y1,polang1,x2,y2,polang2,x3,y3,polang3,x4,y4,polang4)
        write_txt_tmp('FPDB4sim_simfullfpdb_skyview.txt',x0,y0,polang0,x1,y1,polang1,x2,y2,polang2,x3,y3,polang3,x4,y4,polang4)

    if skyview==False:
#        write_txt('simfullfpdb_photonview.txt',x0,y0,polang0,x1,y1,polang1,x2,y2,polang2,x3,y3,polang3,x4,y4,polang4,x5,y5,polang5,x6,y6,polang6)
        write_txt('simfullfpdb_photonview.txt',x0,y0,polang0,x1,y1,polang1,x2,y2,polang2,x3,y3,polang3,x4,y4,polang4)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    q_plot = True
    #q_plot = False
    if q_plot:
        py.plot(x0, y0,'.y')
        py.plot(x0[0],y0[0],'hb')
        py.plot(x0[1],y0[1],'hb')
        py.plot(polvec_x0_t,polvec_y0_t,'r')
        py.plot(polvec_x0_b,polvec_y0_b,'b')
    
        py.plot(x1, y1,'.m')
        py.plot(x1[0],y1[0],'hb')
        py.plot(x1[1],y1[1],'hb')
        py.plot(polvec_x1_t,polvec_y1_t,'r')
        py.plot(polvec_x1_b,polvec_y1_b,'b')
        
        py.plot(x2, y2,'.r')
        py.plot(x2[0],y2[0],'hb')
        py.plot(x2[1],y2[1],'hb')
        py.plot(polvec_x2_t,polvec_y2_t,'r')
        py.plot(polvec_x2_b,polvec_y2_b,'b')    
        
        py.plot(x3, y3,'.g')
        py.plot(x3[0],y3[0],'hb')
        py.plot(x3[1],y3[1],'hb')
        py.plot(polvec_x3_t,polvec_y3_t,'r')
        py.plot(polvec_x3_b,polvec_y3_b,'b')
        
        py.plot(x4, y4,'.k')
        py.plot(x4[0],y4[0],'hb')
        py.plot(x4[1],y4[1],'hb')
        py.plot(polvec_x4_t,polvec_y4_t,'r')
        py.plot(polvec_x4_b,polvec_y4_b,'b')
        
#        py.plot(x5, y5,'.c')
#        py.plot(x5[0],y5[0],'hb')
#        py.plot(x5[1],y5[1],'hb')
#        py.plot(polvec_x5_t,polvec_y5_t,'r')
#        py.plot(polvec_x5_b,polvec_y5_b,'b')
#        
#        py.plot(x6, y6,'.b')
#        py.plot(x6[0],y6[0],'hb')
#        py.plot(x6[1],y6[1],'hb')
#        py.plot(polvec_x6_t,polvec_y6_t,'r')
#        py.plot(polvec_x6_b,polvec_y6_b,'b')
 
        if skyview==False:
            py.xlim([150,-150])
            py.ylim([-150,150])
            py.title('photon view')
        if skyview==True:
            py.xlim([-8,8])
            py.ylim([-8,8])
            py.title('skyview view')
        py.show()


    sys.exit()
###########################################################

    fpdb_list1 = read_FPDB_designed('pixelPosOnly2.txt')
    az1 = fpdb_list1['x_b']
    el1 = fpdb_list1['y_b']
    polang1 = fpdb_list1['polang']
    torb1 = fpdb_list1['torb']
    ind_t1 = np.where('t' == torb1)
    ind_b1 = np.where('b' == torb1)
    ind_t1 = ind_t1[0]; ind_b1 = ind_b1[0]
    polvec_x1, polvec_y1 = gen_polvec(az1,el1,polang1/pi*180.,2)
    
    py.plot(az1,el1,'.')
    py.plot(az1[0],el1[0],'hg')
    py.plot(az1[2],el1[2],'hg')
    for i in ind_t1:
        py.plot(polvec_x1[:,i],polvec_y1[:,i],'b')
    for i in ind_b1:
        py.plot(polvec_x1[:,i],polvec_y1[:,i],'r')
    py.title('')
    py.xlabel('az [degs]')
    py.ylabel('el [degs]')
    # py.savefig(dir+fpdb_fname1+'.png', dpi=DPI)
    
    py.show()


###########################################################
###########################################################
###########################################################

run()
