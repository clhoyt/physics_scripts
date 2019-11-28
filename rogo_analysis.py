# -*- coding: utf-8 -*-
"""
Created on Tue Jun 06 08:22:55 2017

@author: C.L. Hoyt
"""
# Rogowski response

import numpy as np
import matplotlib.pylab as plt
from scipy.constants import inch
plt.close('all')

# constants
eps0 = 8.854e-12
mu0 = 4.0*np.pi*1e-7
rho_cu = 1.68e-8

# RG 223 cable measurements
epsr = 2.25
inner_cond_d = .89e-3
ins_d = 2.95e-3
ins_wall_thk = (ins_d - inner_cond_d)/2.0
wire_d = inch*.0063
wire_r = wire_d/2.0

# cross section of a single turn
A = np.pi*(inner_cond_d)*(inner_cond_d)/4.0 

# number of turns
N =220 # 22

# Rogo diameter
rog_d = 7.122*inch

# Rogo circumferential length
rog_l = np.pi*rog_d

# turns per unit length
n = N/rog_l

# total toroid wire length if unwound
l = np.sqrt((np.pi*ins_d)**2 + (1.0/n)**2)
ll = l*N

dl_ID = 7.122*inch  # dummy load housing inner diameter
a = dl_ID/2.0 - (ins_d+2.0*wire_d)  # for inductance function
b = dl_ID/2.0  # for inductance function

def Rcoil(rho,ll,rog_l,wire_r):
    # rho: resistivity
    # ll: toroid length of wire total
    # rog_l: length of the return wire
    # wire_r: wire radius
    R = rho*(ll+rog_l)/np.pi/wire_r/wire_r
    return R
    

# define coil stray capacitance
def Ccoil(epsr,h,r,l): 
    # h: ground plane to center of conductor distance
    # r: conductor radius
    # l: total length of conductor

    C = 2.0*np.pi*l*eps0*epsr/np.arccosh(h/r)
    return C
    
# define coil inductance
def Lcoil(N,a,b):
    # a: inner major radius
    # b: outer major radius
    # N: number of turns

    L = mu0*N*N/2.0*(a+b-2.0*np.sqrt(a*b))
    return L

# coil mutual coupling     
def Mcoil(N,a,b):
    M = mu0*N/2.0*(a+b-2.0*np.sqrt(a*b))
    return M

def fu(Rc,Rt,Lc,Cc):
    f = (Lc + Rt*Rc*Cc)/(2.0*np.pi*Rt*Lc*Cc)
    return f
    
def fl(Rc,Rt,Lc,Cc):
    f = (Rc+Rt)/(2.0*np.pi*(Lc+Rt*Rc*Cc))
    return f
    
h = ins_wall_thk + wire_r # for capacitance function
r = wire_r # for capacitance function


Lc = Lcoil(N,a,b)
Cc = Ccoil(epsr,h,r,ll)
M = Mcoil(N,a,b)

ww = 1.0/np.sqrt(Lc*Cc)


Rc=Rcoil(rho_cu,ll,rog_l,wire_r)
print("Coil resistance : {}".format(Rc))
Rt=.1
w = np.logspace(1,14,1000)

def Hfunc(Rt,R,L,C,w,M):
    H = -1j*M*w*Rt/(R+Rt-w**2*Rt*L*C+1j*w*(L+R*Rt*C))
    Hmag = np.sqrt(H*H.conjugate())
    return Hmag
    
def Sfunc(Rt,R,L,C,w,M):
    S = w*M*Rt/np.sqrt((R+Rt-w**2*Rt*L*C)**2 + w**2*(L+Rt*R*C)**2)  
    return S    
    
    
def RogCalInt(M,Rt,R,L,C):
    cal = -M*Rt/(L + Rt*R*C)
    return cal

intcal = RogCalInt(M,Rt,Rc,Lc,Cc)
print("Calibration factor: {}".format(intcal))

H = Hfunc(Rt,Rc,Lc,Cc,w,M)
S = Sfunc(Rt,Rc,Lc,Cc,w,M)
flow = fl(Rc,Rt,Lc,Cc)

plt.loglog(w,H,'r',linewidth=2)
plt.loglog(w,S,"c")
plt.axvline(flow,color="red",linestyle="--")
plt.text(2e5,1e-5,"-3dB cutoff",rotation=90)
plt.ylim(1e-7,1e-3)
plt.xlim(1e2,1e10)
#plt.text(1e7,1e-4,"integrating response: 2.2 kA/V \n 10 turns/in. \n 0.1 Ohm resistor to ground")
plt.ylabel("Response (V/A)")
plt.xlabel("Frequency (Hz)")
plt.show()  
    
    
    
    
    
    
    
    
    
    