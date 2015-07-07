# -*- coding: utf-8 -*-
'''
Created on Mon Jul  6 17:43:03 2015

@author: benschneider

System:
SQUID at the end of a Transmission line.
ABCD-Matrix: M1, M2 … ; each represent one element.
For Reference on ABCD Matrix: Microwave Engineering by David M. Pozar p. 185
'''

import numpy as np
from numpy import pi, cos, sin, log
import matplotlib
#matplotlib.use('Qt4Agg')
#matplotlib.use('WX')
matplotlib.use('macosx') # fast on MAC
import matplotlib.pyplot as pl

i = 1.0j
flux0 = 2.07e-15    # Tm^2 ;Flux quanta: flux0 =  h / (2*charging energy)
Z0 = 50.0           # R; Input impedance
Z1 = 1.0           # R; Impedance of cable
Z2 = 50.0           # R; Impedance of Coplanar Waveguide
l1 = 0.1            # m;
v = 2.0e8           # m/s ;approx. velocity in a coaxial 2/3 * speed of light
l2 = 200.0e-6       # m ;adjust length with epsilonr for saphire *3/2
Ic = 0.8e-6         # A ;0.8uA measured, 2.5 uA max
R = 5.0e3           # Ohm
Cap = 20.0e-15      # F
Y4 = 1/0.001           # R ;Wire bonds to GND
xaxis = np.linspace(-1,1,1001)*flux0

#sweep freq
f = 4.0e9           # Hz
b = 2.0*pi*f/v      # 1/m ;b = k = 2pi/wavelength; wavelength = velocity / frequency

#sweep flux
S11_Mag     = np.zeros(xaxis.shape)
S11_Angle   = np.zeros(xaxis.shape)
S12_Mag     = np.zeros(xaxis.shape)
S12_Angle   = np.zeros(xaxis.shape)
ii = 0
for flux in xaxis:
    L = flux0 / (Ic*2.0*pi* np.abs(cos(pi*flux/flux0)))
    Ysq = (1/R + 1/(i*2*pi*f*L +i*1e-90) + i*2*pi*Cap)
    Zsq = 1/Ysq
    leff = L
    l2 = l2 + leff #Squid modulates the phase

    s1 = b*l1
    s2 = b*l2
    #M1 = np.matrix([[cos(s1),i*Z1*sin(s1)],[i*1.0/Z1*sin(s1),cos(s1)]]) # Coaxial Cable with length l1
    #M2 = np.matrix([[cos(s2),i*Z2*sin(s2)],[i*1.0/Z2*sin(s2),cos(s2)]]) # Coplanar Stripline with leght l2
    #M3 = np.matrix([[0,Zsq],[0,1]]) # Perfectly terminated SQUID
    M3 = np.matrix([[1,Zsq],[0,1]]) # Imperfect terminated SQUID
    M4 = np.matrix([[1,0],[Y4,1]]) # Wirebonds to GND
    M = M3*M4
    A = M[0,0]
    B = M[0,1]
    C = M[1,0]
    D = M[1,1]
    S11 = (A+B/Z0-C*Z0-D)/(A+B/Z0+C*Z0+D)
    S12 = 2*(A*D-B*C)/(A+B/Z0+C*Z0+D)

    S11_Mag[ii]     = np.abs(S11)
    S11_Angle[ii]   = np.angle(S11)
    S12_Mag[ii]     = np.abs(S12)
    S12_Angle[ii]   = np.angle(S12)
    ii = ii +1
    #print np.abs(S11)

pl.figure(1)
pl.plot(xaxis/flux,(S11_Mag))
pl.show()