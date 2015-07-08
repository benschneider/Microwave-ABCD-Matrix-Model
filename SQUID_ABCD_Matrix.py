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
from parsers import savemtx, make_header, dim

import matplotlib
matplotlib.use('macosx') # macosx, Qt4Agg, WX
import matplotlib.pyplot as pl

i = 1.0j
flux0 = 2.07e-15    # Tm^2; Flux quanta: flux0 =  h / (2*charging energy)
Z0 = 50.0           # R; Input impedance
Z1 = 20.0           # R; Impedance of cable
Z2 = 50.0           # R; Impedance of Coplanar Waveguide
l1 = 0.44           # m;
l2 = 900.0e-6*2     # m; adjust length with epsilonr for saphire
v = 2.0e8           # m/s; approx. velocity in a coaxial 2/3 * speed of light
Ic = 0.8e-6         # A; 0.8uA measured, 2.5 uA max
R = 5.0e3           # Ohm
Cap = 100.0e-15      # F
Y4 = 1/0.001        # R; Wire bonds to GND


magnet = dim(name = 'Flux (Phi0)',
           start = -1,
           stop = 1,
           pt = 101,
           scale = flux0)
freq = dim(name = 'Frequency (GHz)',
           start = 1,
           stop = 12,
           pt = 1001,
           scale = 1e9)
dim_3 = dim(name = 'Amplitude // Phase',
           start = 0,
           stop = 1,
           pt = 4)
#head1 = make_header(freq, magnet, dim_3, 'S11')
#Mat3d  = np.zeros((dim_3.pt,magnet.pt,freq.pt))

head1 = make_header(magnet, freq, dim_3, 'S11/S12')
Mat3d  = np.zeros((dim_3.pt, freq.pt, magnet.pt))

jj = 0
for f0 in freq.lin:
#for flux in magnet.lin:
    ii = 0
    #for f0 in freq.lin:
    for flux in magnet.lin:
        b = 2.0*pi*f0/v      # b = k = 2pi/wavelength; wavelength = velocity / frequency
        L = flux0 / (Ic*2.0*pi* np.abs(cos(pi*flux/flux0)))
        Ysq = (1.0/R + 1/(i*2*pi*f0*L +i*1e-90) + i*2*pi*Cap)
        Zsq = 1.0/Ysq
        #l2 = l2 + L     #Squid modulates the phase (sin(l) = l for small l) !Wrong!
        s1 = b*l1
        s2 = b*l2
        M1 = np.matrix([[cos(s1),i*Z1*sin(s1)],[i*1.0/Z1*sin(s1),cos(s1)]]) # Coaxial Cable with length l1
        #M3 = np.matrix([[0,1/Ysq],[0,1]]) # Perfectly terminated SQUID
        M2 = np.matrix([[cos(s2),i*Z2*sin(s2)],[i*1.0/Z2*sin(s2),cos(s2)]]) # Coplanar Stripline with leght l2 (including phase modulation of the SQUID)
        M3 = np.matrix([[1,Zsq],[0,1]]) # Non Perfect termination of the SQUID
        M4 = np.matrix([[1,0],[Y4,1]]) # Wirebonds to GND
        M = M1*M2*M3*M4 # connect the elements
        A = M[0,0]
        B = M[0,1]
        C = M[1,0]
        D = M[1,1]
        S11 = (A+B/Z0-C*Z0-D)/(A+B/Z0+C*Z0+D)
        S12 = 2*(A*D-B*C)/(A+B/Z0+C*Z0+D)

        Mat3d[0,jj,ii]     = np.abs(S11)
        Mat3d[1,jj,ii]     = np.angle(S11)
        Mat3d[2,jj,ii]     = np.abs(S12)
        Mat3d[3,jj,ii]     = np.angle(S12)

        ii = ii +1

    Mat3d[1,jj]     = np.unwrap(Mat3d[1,jj])
    Mat3d[3,jj]     = np.unwrap(Mat3d[3,jj])
    jj = jj +1

pl.figure(1)
pl.imshow(Mat3d[0])
pl.show()

pl.figure(2)
pl.imshow(Mat3d[1])
pl.show()

savemtx('resultdata.mtx', Mat3d, header = head1) #mtx file can be opened by spyview
#Link to Spyview: http://nsweb.tn.tudelft.nl/~gsteele/spyview/