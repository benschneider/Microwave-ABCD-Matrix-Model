# -*- coding: utf-8 -*-
'''
Created on Fri Jul 10 15:02:34 2015

@author: benschneider

System:
SQUID at the end of a Transmission line.
ABCD-Matrix: M1, M2 … ; each represent one element.
Ref: 'Microwave Engineering 3rd Edition' by David M. Pozar p. 185
'''
from numpy import pi, cos, log10, abs
import numpy as np
from parsers import dim, get_hdf5data
from ABCD import tline, sres, shunt, handler  # , terminator
# from scipy.io import loadmat, savemat, whosmat #to save and load .mat (matlab)

import matplotlib
matplotlib.use('Qt4Agg')  # macosx, Qt4Agg, WX
# matplotlib.use('macosx')  # macosx, Qt4Agg, WX
import matplotlib.pyplot as plt

flux0 = 2.07e-15    # Tm^2; Flux quanta: flux0 =  h / (2*charging energy)

squid = dim(name='Squid',
            start=-1.5,
            stop=0.5,
            pt=501,
            scale=flux0)
elem = handler(name='mag/phase',
               start=0,
               stop=10,
               pt=1)  # 8 pts for S 4x2 values
# squid.Ic=1.7e-6       # A; Ic ~ 0.85uA measured, 2.5 uA max
# squid.R = 2.3e3          # Ohm
# squid.Cap = 450e-15     # 450.0e-15     # F
squid.flux0 = 2.07e-15  # Tm^2; Flux quanta: flux0 =  h / (2*charging energy)

elem.Z0 = 50            # R; Input impedance
elem.Z1 = 50            # R; Impedance of transmission piece 1
elem.Z2 = 50            # R; Impedance of Coplanar Waveguide
elem.Z3 = 50
elem.L1 = 0.44
elem.L2 = 900.0e-6
elem.L3 = 0.01
elem.Z4 = 0.1           # Ohm; Wire bonds conductance to GND (-45dB isolation)
# m/s; approx. velocity in a coaxial 2/3 * speed of light
elem.v = 2.0e8

# load hdf5 data
# filename = 'S1_203.hdf5'
filename = 'S1_945_S11_4p1_BPF7.hdf5'
# filename = 'S1_945_S11_4p8_BPF7.hdf5'
measdata = get_hdf5data(filename)


def get_sMatrix(b, elem, Zsq):
    SM = np.zeros((len(Zsq), 2, 2))*1j  # complex matrix
    M1 = (tline(elem.Z1, b, elem.L1) *
          tline(elem.Z2, b, elem.L2) *
          tline(elem.Z3, b, elem.L3))  # transmission lines
    for ii, Zsq1 in enumerate(Zsq):
        M2 = sres(Zsq1)*shunt(elem.Z4)
        M4 = M1*M2
        SM[ii] = elem.get_SM(M4)  # complex S-Matrix shape [2 2 fluxlength]
    return SM


def get_Zsq(f0, squid):
    flux = squid.lin
    L = flux0 / (squid.Ic * 2.0 * pi * abs(cos(pi*flux / squid.flux0)))
    Ysq = (1.0/squid.R
           + 1.0/(1j*2.0*pi*f0*L + 1j*1e-90)
           + 1j*2.0*pi*f0*squid.Cap)
    return (1.0/Ysq)


def get_SMresponse(f0, squid, elem):
    b = 2.0*pi*f0/2.0e-8
    Zsq = get_Zsq(f0, squid)
    return get_sMatrix(b, elem, Zsq)


plt.ion()


def plotfig2(SMat, measdata):
    ydat = measdata.D1complex[:, measdata.findex]
    S11 = SMat[:, 0, 0]
    # S12 = SMat[:,1,0]
    xaxis = squid.lin/flux0
    xaxis2 = np.linspace(-1+measdata.XPOS,
                         1+measdata.XPOS,
                         len(ydat))*measdata.XSC
    fig2 = plt.figure(2)
    g1 = fig2.add_subplot(2, 1, 1)
    g1.plot(xaxis, abs(S11) * 10**(measdata.ATT/20.0))
    g1.hold(True)
    g1.plot(xaxis2, abs(ydat))
    # g1.set_ylim([0.9,1.0])
    g1.hold(False)
    g2 = fig2.add_subplot(2, 1, 2)  # , sharex=g1)
    g2.plot(xaxis, np.unwrap(np.angle(S11))*180/pi + measdata.PHI)
    g2.hold(True)
    g2.plot(xaxis2, np.unwrap(np.angle(ydat), discont=pi/2) * 180/pi)
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


def plotfig4(SMat, measdata):
    ydat = measdata.D1complex[:, measdata.findex]
    xaxis2 = np.linspace(-1+measdata.XPOS,
                         1+measdata.XPOS,
                         len(ydat))*measdata.XSC
    fig4 = plt.figure(4)
    g1 = fig4.add_subplot(2, 1, 1)
    g1.plot(xaxis2, abs(ydat))
    g1.hold(False)
    g2 = fig4.add_subplot(2, 1, 2)  # , sharex=g1)
    g2.plot(xaxis2, np.unwrap(np.angle(ydat), discont=pi/2) * 180/pi)
    g2.hold(False)
    fig4.show()


def fitcurve(val):
    #  take current data to use of a 2d fit
    print squid.Ic

execfile('interface.py')
