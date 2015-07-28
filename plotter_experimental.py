# -*- coding: utf-8 -*-
'''
Created on Fri Jul 10 15:02:34 2015

@author: benschneider

System:
SQUID at the end of a Transmission line.
ABCD-Matrix: M1, M2 … ; each represent one element.
Ref: 'Microwave Engineering 3rd Edition' by David M. Pozar p. 185
'''

<<<<<<< HEAD
from numpy import pi, cos, log10, abs
import numpy as np
from parsers import dim, get_hdf5data
=======
from numpy import pi, cos, abs, zeros, angle, unwrap, log #, sin, log
from parsers import dim #make_header, savemtx
>>>>>>> origin/master
from ABCD import tline, sres, shunt, handler, terminator
#from scipy.io import loadmat, savemat, whosmat #to save and load .mat (matlab)

import matplotlib
matplotlib.use('Qt4Agg') # macosx, Qt4Agg, WX
import matplotlib.pyplot as plt

flux0 = 2.07e-15    # Tm^2; Flux quanta: flux0 =  h / (2*charging energy)

squid = dim(name = 'Squid',
            start = -1,
            stop = 1,
            pt = 216,
            scale = flux0)
elem = handler(name = 'mag/phase',
               start = 0,
               stop = 10,
               pt = 1) #8 pts for S 4x2 values
#squid.Ic = 1.7e-6       # A; Ic ~ 0.85uA measured, 2.5 uA max
#squid.R = 2.3e3          # Ohm
#squid.Cap = 450e-15     # 450.0e-15     # F
squid.flux0 = 2.07e-15  # Tm^2; Flux quanta: flux0 =  h / (2*charging energy)

elem.Z0 = 50            # R; Input impedance
elem.Z1 = 50            # R; Impedance of transmission piece 1
elem.Z2 = 50            # R; Impedance of Coplanar Waveguide
elem.Z3 = 50
elem.L1 = 0.44
elem.L2 = 900.0e-6
elem.L3 = 0.01
elem.Z4 = 0.1           # Ohm; Wire bonds conductance to GND (-45dB isolation)
elem.v = 2.0e8          # m/s; approx. velocity in a coaxial 2/3 * speed of light

#load hdf5 data
filename = 'S1_203.hdf5'
measdata = get_hdf5data(filename)

def get_sMatrix(b,elem,Zsq):
    SM = np.zeros( (len(Zsq), 2, 2) )*1j #complex matrix
    M1 = (tline(elem.Z1,b,elem.L1)*
            tline(elem.Z2,b,elem.L2)*
            tline(elem.Z3,b,elem.L3)) #transmission lines
    for ii, Zsq1 in enumerate(Zsq):
        M2 = sres(Zsq1)*shunt(elem.Z4)
        M4 = M1*M2
        SM[ii] = elem.get_SM(M4) #complex S-Matrix shape [2 2 fluxlength]
    return SM

def get_Zsq(f0,squid):
    flux = squid.lin
    L = flux0 / (squid.Ic*2.0*pi* abs(cos(pi*flux/squid.flux0)))
    Ysq = (1.0/squid.R + 1.0/(1j*2.0*pi*f0*L + 1j*1e-90) + 1j*2.0*pi*f0*squid.Cap)
    return (1.0/Ysq)

def get_SMresponse(f0,squid,elem):
    b = 2.0*pi*f0/2.0e-8
    Zsq = get_Zsq(f0, squid)
    return get_sMatrix(b,elem,Zsq)

def find_nearest(someArray,value):
    idx = abs(someArray-value).argmin()
    return idx

plt.ion()

def plotfig2(SMat,measdata):
    ydat = measdata.D1complex[:,measdata.findex]
    S11 = SMat[:,0,0]
    #S12 = SMat[:,1,0]
    xaxis = squid.lin/flux0
    xaxis2 = np.linspace(-1+measdata.XPOS,1+measdata.XPOS,len(ydat))*measdata.XSC
    fig2 = plt.figure(2)
<<<<<<< HEAD
    g1 = fig2.add_subplot(2, 1, 1)
    g1.plot(xaxis, 20*log10(abs(S11))+ measdata.ATT )
    g1.hold(True)
    g1.plot(xaxis2, 20*log10(abs(ydat)))
=======
    g1 = fig2.add_subplot(2, 2, 1)
    g1.plot(xaxis, (abs(S11)))
>>>>>>> origin/master
    #g1.set_ylim([0.9,1.0])
    g1.hold(False)
    g2 = fig2.add_subplot(2, 1, 2) #, sharex=g1)
    g2.plot(xaxis, np.unwrap(np.angle(S11))*180/pi + measdata.PHI)
    g2.hold(True)
    g2.plot(xaxis2, np.unwrap(np.angle(ydat), discont = pi/2) *180/pi)
    g2.hold(False)
    '''
    g3 = fig2.add_subplot(2, 2, 2, sharex=g1)
    g3.plot(xaxis, abs(S12))
    g3.hold(False)
    g4 = fig2.add_subplot(2, 2, 4, sharex=g1)
    g4.plot(xaxis, unwrap(angle(S12))*180/pi)
    g4.hold(False)
    '''
    fig2.show()

def plotfig4(SMat,measdata):
    ydat = measdata.D1complex[:,measdata.findex]
    xaxis2 = np.linspace(-1+measdata.XPOS,1+measdata.XPOS,len(ydat))*measdata.XSC
    fig4 = plt.figure(4)
    g1 = fig4.add_subplot(2, 1, 1)
    g1.plot(xaxis2, 20*log10(abs(ydat)))
    g1.hold(False)
    g2 = fig4.add_subplot(2, 1, 2) #, sharex=g1)
    g2.plot(xaxis2, np.unwrap(np.angle(ydat), discont = pi/2) *180/pi)
    g2.hold(False)
    fig4.show()



fig3 = plt.figure(3)
axcolor = 'lightgoldenrodyellow'
axfreq = plt.axes([0.25, 0.1, 0.50, 0.03], axisbg=axcolor)
axIc  = plt.axes([0.25, 0.15, 0.50, 0.03], axisbg=axcolor)
axCap  = plt.axes([0.25, 0.20, 0.50, 0.03], axisbg=axcolor)
axRsq  = plt.axes([0.25, 0.25, 0.50, 0.03], axisbg=axcolor)
axZ1  = plt.axes([0.25, 0.30, 0.50, 0.03], axisbg=axcolor)
axL1  = plt.axes([0.25, 0.35, 0.50, 0.03], axisbg=axcolor)
axZ2  = plt.axes([0.25, 0.40, 0.50, 0.03], axisbg=axcolor)
axL2  = plt.axes([0.25, 0.45, 0.50, 0.03], axisbg=axcolor)
axZ3  = plt.axes([0.25, 0.50, 0.50, 0.03], axisbg=axcolor)
axL3  = plt.axes([0.25, 0.55, 0.50, 0.03], axisbg=axcolor)
axZ4  = plt.axes([0.25, 0.60, 0.50, 0.03], axisbg=axcolor)
axXPOS  = plt.axes([0.25, 0.70, 0.50, 0.03], axisbg=axcolor)
axXSC  = plt.axes([0.25, 0.75, 0.50, 0.03], axisbg=axcolor)
axATT  = plt.axes([0.25, 0.80, 0.50, 0.03], axisbg=axcolor)
axPHI  = plt.axes([0.25, 0.85, 0.50, 0.03], axisbg=axcolor)

sFreq = plt.Slider(axfreq, 'Freq (GHz)', measdata.freq[0]/1e9, measdata.freq[-1]/1e9, valinit=5.3)
sIc = plt.Slider(axIc, 'Ic (uA)', 0.1, 10.0, valinit=3.3)
sCap = plt.Slider(axCap, 'Cap (fF)', 0, 500.0, valinit=40)
sRsq = plt.Slider(axRsq, 'Rsq (kOhm)', 0.01, 10.0, valinit=5)
sZ1 = plt.Slider(axZ1, 'Z1 (Ohm)', 0.0, 900.0, valinit=50)
sZ2 = plt.Slider(axZ2, 'Z2 (Ohm)', 0.0, 400.0, valinit=50)
sZ3 = plt.Slider(axZ3, 'Z3 (Ohm)', 0.0, 400.0, valinit=50)
sL1 = plt.Slider(axL1, 'L1 (m)', 0.0, 0.1, valinit=0.01)
sL2 = plt.Slider(axL2, 'L2 (m)', 0.0, 1.0, valinit=0.3)
sL3 = plt.Slider(axL3, 'L3 (mm)', 0.0, 20.0, valinit=0.9)
sZ4 = plt.Slider(axZ4, 'W.b. -> GND (Ohm)', 0.0001, 1.0, valinit=0.1)
sXPOS = plt.Slider(axXPOS, 'x-pos', -1.0, 1.0, valinit=0)
sXSC = plt.Slider(axXSC, 'x-scale', 0.5, 1.5, valinit=0.9)
sATT = plt.Slider(axATT, 'Attenuation dBm', -70, 0.0, valinit=-50)
sPHI = plt.Slider(axPHI, 'Phase offset', -360, 360.0, valinit= 0)

fig3.show()

def update(val):
    f0 = sFreq.val*1e9
    measdata.findex = find_nearest(measdata.freq,f0)
    squid.Ic = sIc.val*1e-6
    squid.Cap = sCap.val*1e-15
    squid.R = sRsq.val*1e3
    elem.Z1 = sZ1.val
    elem.L1 = sL1.val
    elem.Z2 = sZ2.val
    elem.L2 = sL2.val
    elem.Z3 = sZ3.val
    elem.L3 = sL3.val*1e-3
    elem.Z4 = sZ4.val
    measdata.XSC = sXSC.val
    measdata.XPOS = sXPOS.val
    measdata.ATT = sATT.val
    measdata.PHI = sPHI.val
    SMat = get_SMresponse(f0,squid,elem)
    plotfig2(SMat,measdata)
    plotfig4(SMat,measdata)


sIc.on_changed(update)
sCap.on_changed(update)
sRsq.on_changed(update)
sZ1.on_changed(update)
sZ2.on_changed(update)
sZ3.on_changed(update)
sL1.on_changed(update)
sL2.on_changed(update)
sL3.on_changed(update)
sZ4.on_changed(update)
sXPOS.on_changed(update)
sXSC.on_changed(update)
sATT.on_changed(update)
sPHI.on_changed(update)
sFreq.on_changed(update)

updatetax = plt.axes([0.1, 0.025, 0.1, 0.04])
button2 = plt.Button(updatetax, 'Update', color=axcolor, hovercolor='0.975')
button2.on_clicked(update)

update(0)
