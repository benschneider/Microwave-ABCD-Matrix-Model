# -*- coding: utf-8 -*-
'''
Created on Mon Jul  6 17:43:03 2015

@author: benschneider

System:
SQUID at the end of a Transmission line.
ABCD-Matrix: M1, M2 … ; each represent one element.
Ref: 'Microwave Engineering 3rd Edition' by David M. Pozar p. 185
'''
import numpy as np
from numpy import pi, cos, sin, log
from parsers import savemtx, make_header, dim
from ABCD import *


import matplotlib
matplotlib.use('macosx') # macosx, Qt4Agg, WX
import matplotlib.pyplot as pl

i = 1.0j
flux0 = 2.07e-15    # Tm^2; Flux quanta: flux0 =  h / (2*charging energy)
Z0 = 50.0           # R; Input impedance
Z1 = 30.0           # R; Impedance of transmission piece 1
Z2 = 50.0           # R; Impedance of Coplanar Waveguide
l1 = 0.44           # m;
l2 = 900.0e-6       # m; adjust length with epsilonr for saphire
v = 2.0e8           # m/s; approx. velocity in a coaxial 2/3 * speed of light
Ic = 1.7e-6         # A; Ic ~ 0.85uA measured, 2.5 uA max
R = 2.3e3           # Ohm
Cap = 450e-15       # 450.0e-15     # F
Y4 = 1/0.1          # 1/Ohm; Wire bonds conductance to GND (-45dB isolation)

magnet = dim(name = 'Flux (Phi0)',
           start = -1,
           stop = 1,
           pt = 1001,
           scale = flux0)
freq = dim(name = 'Frequency (GHz)',
           start = 4,
           stop = 8,
           pt = 101,
           scale = 1e9)
dim_3 = handler(name = 'mag/phase',
           start = 0,
           stop = 10,
           pt = 11) #8 pts for S 4x2 values
dim_3._Z0 = 50

head1 = make_header(magnet, freq, dim_3, 'S11 S12 S21 S22 Z L')
dim_3.prepare_data_save(magnet, freq, dim_3)


jj = 0
for f0 in freq.lin:
    ii = 0
    for flux in magnet.lin:
        # b = k = 2pi/wavelength; wavelength = velocity / frequency
        b = 2.0*pi*f0/v
        L = flux0 / (Ic*2.0*pi* np.abs(cos(pi*flux/flux0)))
        Ysq = (1.0/R + 1.0/(i*2.0*pi*f0*L +i*1e-90) + i*2.0*pi*f0*Cap)
        Zsq = 1.0/Ysq

        ABCD = tline(70,b,0.01)*tline(50,b,0.3)*tline(10,b,900e-6)*sres(Zsq)*shunt(0.1)

        #record stuff into dim_3._SMat
        dim_3.record_SM(ABCD,jj,ii)
        dim_3.record_ZL(Zsq,L, jj,ii)

        ii = ii +1
    dim_3.unwrap_SM(jj)
    dim_3._SMat[9,jj] = np.unwrap(dim_3._SMat[9,jj])
    jj = jj +1

pl.figure(1)
pl.imshow(dim_3._SMat[0])
pl.show()

pl.figure(2)
pl.imshow(dim_3._SMat[1])
pl.show()

savemtx('resultdata3.mtx', dim_3._SMat, header = head1) #mtx file can be opened by spyview
#Link to Spyview: http://nsweb.tn.tudelft.nl/~gsteele/spyview/