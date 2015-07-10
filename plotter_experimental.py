# -*- coding: utf-8 -*-
'''
Created on Fri Jul 10 15:02:34 2015

@author: benschneider

System:
SQUID at the end of a Transmission line.
ABCD-Matrix: M1, M2 … ; each represent one element.
Ref: 'Microwave Engineering 3rd Edition' by David M. Pozar p. 185
'''

from numpy import pi, cos, abs, unwrap #, sin, log
from parsers import savemtx, make_header, dim
from ABCD import handler, tline, sres, shunt


import matplotlib
matplotlib.use('macosx') # macosx, Qt4Agg, WX
import matplotlib.pyplot as plt

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


for jj, f0 in enumerate(freq.lin):
    for ii, flux in enumerate(magnet.lin):
        # b = k = 2pi/wavelength; wavelength = velocity / frequency
        b = 2.0*pi*f0/v
        L = flux0 / (Ic*2.0*pi* abs(cos(pi*flux/flux0)))
        Ysq = (1.0/R + 1.0/(i*2.0*pi*f0*L +i*1e-90) + i*2.0*pi*f0*Cap)
        Zsq = 1.0/Ysq

        ABCD_Matrix = tline(70,b,0.01)*tline(50,b,0.3)*tline(10,b,900e-6)*sres(Zsq)*shunt(0.1)

        #record stuff into dim_3._SMat
        dim_3.record_SM(ABCD_Matrix,jj,ii)
        dim_3.record_ZL(Zsq,L, jj,ii)

    dim_3.unwrap_SM(jj)
    dim_3._SMat[9,jj] = unwrap(dim_3._SMat[9,jj])

plt.figure(1)
plt.subplot(2, 1, 1)
plt.imshow(dim_3._SMat[0], aspect = 'auto',cmap=plt.get_cmap('seismic'))
plt.subplot(2, 1, 2)
plt.imshow(dim_3._SMat[1], aspect = 'auto',cmap=plt.get_cmap('seismic'))
plt.show()


def plotfig2(f):
    plt.figure(2)
    plt.subplot(2, 1, 1)
    plt.plot(dim_3._SMat[0][f])
    plt.hold(False)
    plt.subplot(2, 1, 2)
    plt.plot(dim_3._SMat[1][f])
    plt.hold(False)
    plt.show()


plt.figure(3)
axcolor = 'lightgoldenrodyellow'
axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
axamp  = plt.axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
sfreq = plt.Slider(axfreq, 'Freq', 0, 100.0, valinit=1)
samp = plt.Slider(axamp, 'Amp', 0.1, 10.0, valinit=1)

def update(val):
    #amp = samp.val
    freq = int(sfreq.val)
    plotfig2(freq)
    #plt.hold(False)
sfreq.on_changed(update)
samp.on_changed(update)


#savemtx('resultdata3.mtx', dim_3._SMat, header = head1) #mtx file can be opened by spyview
#Link to Spyview: http://nsweb.tn.tudelft.nl/~gsteele/spyview/