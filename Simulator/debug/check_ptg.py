import numpy as np
import pylab as py
import os

pi = np.pi
radeg = (180./pi)

out = np.load('tmp.npz')
i_pix = out['i_pix']
pix_list = out['pix_list']
ptg_package = out['ptg_package']
fpdb_list_simgen = out['fpdb_list_simgen']
fpdb_list_mmin = out['fpdb_list_mmin']
top_ptg_in = out['top_ptg_in']
bot_ptg_in = out['bot_ptg_in']
top_ptg_out = out['top_ptg_out']
bot_ptg_out = out['bot_ptg_out']

num1=0;num2=100000
bs_glat = (ptg_package.item())['glat'][num1:num2]*radeg
bs_glon = (ptg_package.item())['glon'][num1:num2]*radeg
bs_pa = (ptg_package.item())['pa'][num1:num2]*radeg

t_glat = (top_ptg_in.item())['glat'][num1:num2]*radeg
t_glon = (top_ptg_in.item())['glon'][num1:num2]*radeg
t_psi = (top_ptg_in.item())['psi'][num1:num2]*radeg

b_glat = (bot_ptg_in.item())['glat'][num1:num2]*radeg
b_glon = (bot_ptg_in.item())['glon'][num1:num2]*radeg
b_psi = (bot_ptg_in.item())['psi'][num1:num2]*radeg

py.figure(figsize=(15,10))
py.subplot(241)
py.plot(bs_glat,'.')
py.title('bs_glat')

py.subplot(242)
py.plot(bs_glon,'.')
py.title('bs_glon')

py.subplot(243)
py.plot(bs_glon,bs_glat,'.')
py.title('bs_glon/bs_glat')

py.subplot(244)
py.plot(bs_pa,'.')
py.title('bs_pa')

py.subplot(245)
py.plot(t_glat,'x')
py.plot(b_glat,'.')
py.title('t_glat')

py.subplot(246)
py.plot(t_glon,'x')
py.plot(b_glon,'.')
py.title('t_glon')

py.subplot(247)
py.plot(t_glon,t_glat,'x')
py.plot(b_glon,b_glat,'.')
py.title('t_glon/t_glat')

py.subplot(248)
py.plot(t_psi,'.')
py.plot(b_psi,'.')
py.title('t and b_psi')

py.savefig('tmp.png')
os.system('display tmp.png &')
