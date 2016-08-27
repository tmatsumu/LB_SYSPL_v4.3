from numpy import *
from pylab import *
from toollib_common import *
from constant import *
import sys

# Plot circle or radius 3
def GaussBeam(x,y,x0,y0,fwhm):
    sigma = fwhm/(2.*sqrt(2.*log(2.)))
    sigma = sigma*0.7
    out = exp((-(x-x0)**2-(y-y0)**2)/(2.*sigma**2.))
    ind = where( (x-x0)**2+(y-y0)**2 > (0.5*fwhm*1.1)**2 )
    out[ind] = 0.
    return out

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c = c_m() # in mm
pi = pi()
r2d = r2d()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

fp_inside = float(sys.argv[1])
D_lens1 = float(sys.argv[2])
D_sep1 = float(sys.argv[3])
fp_outside_x = float(sys.argv[4])
fp_outside_y = float(sys.argv[5])
D_lens2 = float(sys.argv[6])
D_sep2 = float(sys.argv[7])
f_eff = float(sys.argv[8])
freq = float(sys.argv[9])
D_aperture = float(sys.argv[10])*1e-2 # cm->m
plot_type = sys.argv[11]
in_out = sys.argv[12]

print '++++++++++++++++++++++++++++++++++++'
print 'plot type: ', plot_type
print 'Frequency: ', freq*1e-9, 'GHz'
if 'FP' in plot_type:
    plate_scale = 1.
    wire = 1.2
if 'FOV' in plot_type:
    plate_scale = 1./f_eff*180./pi
    wire = 1.
    fwhm_beam1 =2.* arccos( 1.-(c/freq)**2.*(2./D_lens1)**2./(2.*pi**2) ) *r2d
    fwhm_beam2 =2.* arccos( 1.-(c/freq)**2.*(2./D_lens2)**2./(2.*pi**2) ) *r2d
    print 'FWHM beam 1 [arcmin]: ', fwhm_beam1*60.
    print 'FWHM beam 2 [arcmin]: ', fwhm_beam2*60.
    D_lens1 =2.* arccos( 1.-(c/freq)**2.*(2./D_aperture)**2./(2.*pi**2) ) *r2d
    D_lens2 = D_lens1
    print 'FWHM [arcmin]: ', D_lens1*60.

D_sep1 = D_sep1*plate_scale
D_sep2 = D_sep2*plate_scale
fp_outside_x = fp_outside_x*plate_scale
fp_outside_y = fp_outside_y*plate_scale
Dfp1 = fp_inside*plate_scale

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# make these smaller to increase the resolution
dx, dy = 0.05, 0.05

x = arange(-fp_outside_x/2.*1.2, fp_outside_x/2.*1.2, dx)
y = arange(-fp_outside_y/2.*1.2, fp_outside_y/2.*1.2, dy)
X,Y = meshgrid(x, y)
Zout = GaussBeam(X,Y,0.,0.,1.,) * 0.

an = linspace(0,2*pi,100)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

R_lens1 = D_lens1/2.
R_sep1 = D_sep1/2.

R_sep2 = D_sep2/2.
Dfp2_x = fp_outside_x
Dfp2_y = fp_outside_y
Dfp2 = Dfp2_x

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
f_txt = open("test.txt", "w")
if ('in' in in_out) | ('both' in in_out):
    y0_tmp=arange(Dfp1/(2.*R_sep1)*1.5)*sqrt(3.)*R_sep1
    y0_tmpr=y0_tmp[1:Dfp1/(2.*R_sep1)*1.5-1]
    y0=concatenate((y0_tmp,-y0_tmpr[::-1]))
    nby=len(y0)

    count1=0
    for iy in range(0,nby):
        y0in=y0[iy]
        amari1 = y0in/(sqrt(3.)*R_sep1*2.) - int(y0in/(sqrt(3.)*R_sep1*2.))

        if (abs(amari1) < 1e-10):
            for ix in range(-nby,nby):
                x0in = ix*2.*R_sep1
                if (sqrt(x0in**2.+y0in**2.) < Dfp1/2.):
                    if 'FP' in plot_type: 
                        plot( 0.5*D_lens1*cos(an)+x0in, 0.5*D_lens1*sin(an)+y0in, 'r')
                        f_txt.write('%d %f %f \n' % (iy, x0in, y0in))
                    if 'FOV' in plot_type: Zout += GaussBeam(X, Y, x0in, y0in, D_lens1)
                    count1+=1
        if (abs(amari1) > 0.48):
            for ix in range(-nby,nby):
                x0in = (ix*2.*R_sep1+R_sep1)
                if (sqrt(x0in**2.+y0in**2.) < Dfp1/2.):
                    if 'FP' in plot_type: 
                        plot( 0.5*D_lens1*cos(an)+x0in, 0.5*D_lens1*sin(an)+y0in, 'r')
                        f_txt.write('%d %f %f \n' % (iy, x0in, y0in))
                    if 'FOV' in plot_type: Zout += GaussBeam(X, Y, x0in, y0in, D_lens1) 
                    count1+=1

    print 'Npix inside=', count1

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if ('out' in in_out) | ('both' in in_out):
    y0_tmp=arange(Dfp2/(2.*R_sep2)*1.5)*sqrt(3.)*R_sep2
    y0_tmpr=y0_tmp[1:Dfp2/(2.*R_sep2)*1.5-1]
    y0=concatenate((y0_tmp,-y0_tmpr[::-1]))
    nby=len(y0)

    row_sign = -1
    count2=0
    count_y=0
    for iy in range(0,nby):
        y0in=y0[iy]

        amari2 = y0in/(sqrt(3.)*R_sep2*2.) - int(y0in/(sqrt(3.)*R_sep2*2.))
        if ((abs(amari2) > 0.1) and (abs(amari2) < 0.9)): amari2=0.5
        if ((abs(amari2) < 0.1) | (abs(amari2) > 0.9)): amari2=0.

        # row of amari2 = 0
        if (abs(amari2) < 1e-5):
            row_sign *= -1

            for ix in range(-nby*2,nby*2):
                x0in = ix*2.*R_sep2

                # draw a Gaussian beam within the outside ellipse, leave the center hole empty
                if ((sqrt((x0in/(Dfp2_x/2.))**2.+(y0in/(Dfp2_y/2.))**2.) < 1.) and (sqrt((x0in/(Dfp1*1.2/2.))**2.+(y0in/(Dfp1*1.2/2.))**2.) > 1.)):
                    if 'FP' in plot_type:
                        plot( 0.5*D_lens2*cos(an)+x0in, 0.5*D_lens2*sin(an)+y0in, 'b')
                        f_txt.write('%d %f %f \n' % (iy, x0in, y0in))
                        count2+=1
                    if 'FOV' in plot_type:
                        if (freq > 90e9):
                            Zout += GaussBeam(X, Y, x0in, y0in, D_lens2)
                            count2+=1
                        if (freq < 90e9):
#                            print 'x, y', ix/2.-int(ix/2), abs(count_y/2.-int(count_y/2))
                            # place the beam in every other spot
                            if (ix/2. - int(ix/2)) == 0:
                                if row_sign == 1: Zout += GaussBeam(X, Y, x0in, y0in, D_lens2)
                                if ((row_sign == -1) and ((sqrt( ((x0in+R_sep2*2)/(Dfp2_x/2.))**2.+(y0in/(Dfp2_y/2.))**2. ) < 1.) and (sqrt( ((x0in+R_sep2*2)/(Dfp1*1.2/2.))**2.+(y0in/(Dfp1*1.2/2.))**2.) > 1.))):
                                    Zout += GaussBeam(X, Y, x0in+R_sep2*2, y0in, D_lens2)
#                                if abs(count_y/2.-int(count_y/2)) < 0.1: Zout += GaussBeam(X, Y, x0in, y0in, D_lens2)
#                                if abs(count_y/2.-int(count_y/2)) > 0.48: Zout += GaussBeam(X, Y, x0in+R_sep2*2, y0in, D_lens2)
                                count2+=1
                                count_y+=1

        if (abs(amari2) > 0.48):
            for ix in range(-nby*2,nby*2):
                x0in = (ix*2.*R_sep2+R_sep2)
                #            if ((sqrt(x0in**2.+y0in**2.) < Dfp2/2.) and (sqrt(x0in**2.+y0in**2.)) > Dfp1/2.):
                if ((sqrt((x0in/(Dfp2_x/2.))**2.+(y0in/(Dfp2_y/2.))**2.) < 1.) and (sqrt((x0in/(Dfp1*1.2/2.))**2.+(y0in/(Dfp1*1.2/2.))**2.) > 1.)):
                    if 'FP' in plot_type:
                        plot( 0.5*D_lens2*cos(an)+x0in, 0.5*D_lens2*sin(an)+y0in, 'b')
                        f_txt.write('%d %f %f \n' % (iy, x0in, y0in))
                        count2+=1
                    if 'FOV' in plot_type:
                        if (freq > 90e9):
                            Zout += GaussBeam(X, Y, x0in, y0in, D_lens2)
                            count2+=1

    print 'Npix outside=', count2

f_txt.close()

#title('not equal, looks like ellipse',fontsize=10)
#xlim(-10,10)
#ylim(-10,10)

axis('equal')
if 'FOV' in plot_type:
    print 'plate scale', plate_scale
#    pcolor(X, Y*plate_scale, Zout)
#    axis([-fp_outside_x/2.*plate_scale,fp_outside_x/2.*plate_scale,-fp_outside_y/2.*plate_scale,fp_outside_y/2.*plate_scale])
#    ind = where(Zout > 1.2)
#    Zout[ind] = 0.
    pcolor(X, Y, Zout)
    axis([-fp_outside_x/2.,fp_outside_x/2.,-fp_outside_y/2.,fp_outside_y/2.])
    colorbar()

show()
