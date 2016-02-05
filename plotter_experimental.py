# -*- coding: utf-8 -*-
'''
Created on Fri Jul 10 15:02:34 2015
@author: benschneider
System:
SQUID at the end of a Transmission line.
ABCD-Matrix: M1, M2 … ; each represent one element.
Ref: 'Microwave Engineering 3rd Edition' by David M. Pozar p. 185
'''
from numpy import pi, cos, abs, log10
import numpy as np
from time import time
from parsers import dim, get_hdf5data
from ABCD import tline, sres, shunt, handler  # , terminator
# from scipy.optimize import curve_fit  # , leastsq
# from scipy.io import loadmat, savemat, whosmat #to save and load .mat (matlab)
from lmfit import minimize, Parameters, Parameter, report_fit
import matplotlib
matplotlib.use('Qt4Agg')  # macosx, Qt4Agg, WX
import matplotlib.pyplot as plt
# matplotlib.use('macosx')

plt.ion()

# load hdf5 data
# filename = 'S1_203.hdf5'
filename = 'S1_945_S11_4p1_BPF7.hdf5'
# filename = 'S1_945_S11_4p8_BPF7.hdf5'
measdata = get_hdf5data(filename)

flux0 = 2.07e-15    # Tm^2; Flux quanta: flux0 =  h / (2*charging energy)

squid = dim(name='Squid',
            start=-1.5,
            stop=0.5,
            pt=201,
            scale=flux0)
squid.matchX = False
elem = handler(name='mag/phase',
               start=0,
               stop=10,
               pt=1)  # 8 pts for S 4x2 values
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
squid.Lbw = 1e-90
squid.Lloop = 1e-90

elem.update = True      # Set update to true


def get_sMatrix(b, elem, Zsq):
    '''
    This is where the assembly of the model happens using ABCD Matrixes
    SM represents the SQUID response
    M1 respresents 3 sections of Transmission lines (L1,Z1)(L2,Z2)(L3,Z3)
    Zsq is the impedance of the SQUID (flux) (vector)
    Assembly is done for each flux point and change in Zsq
    '''
    SM = np.zeros((len(Zsq), 2, 2)) * 1j  # complex matrix
    # M1 = (tline(elem.Z1, b, elem.L1) *
    #       tline(elem.Z2, b, elem.L2) *
    #       tline(elem.Z3, b, elem.L3))  # transmission lines
    M1 = tline(elem.Z1, b, elem.L1)
    for ii, Zsq1 in enumerate(Zsq):
        M2 = sres(Zsq1) * shunt(elem.Z4)
        M4 = M1 * M2
        SM[ii] = elem.get_SM(M4)  # complex S-Matrix shape [2 2 fluxlength]
    return SM


def get_Zsq(f0, squid):
    '''
    Calculates the Impedance of the SQUID.
    First calculates the inductance L  (flux)
    Then calulcates the Impedance (Cap, Freq, L)
    '''
    flux = squid.lin
    L = flux0/(squid.Ic*2.0*pi*abs(cos(pi*flux/squid.flux0)))
    Ysq = (1.0/squid.R + 1j*2.0*pi*f0*squid.Cap - 1j/(2.0*pi*f0*L+squid.Lbw))
    return (1.0/Ysq)


def get_SMresponse(f0, squid, elem):
    '''
    b is a wavenumber used for the Transmission line sections
    2e8 is ~ speed of light in a transmission line
    Zsq vector is the Squid impedance (flux)
    get_sMatrix assembles the individual Transmission sections
    '''
    b = 2.0 * pi * f0 / 2.0e-8
    Zsq = get_Zsq(f0, squid)
    return get_sMatrix(b, elem, Zsq)


def getModelData(squid, elem, measdata):
    '''
    Get S11 including Attenuation and Phase
    is using the full ABCD Model
    Phase and Attenuation is added (outside the ABCD model)
    '''
    f0 = sFreq.val * 1e9
    SMat = get_SMresponse(f0, squid, elem)
    S11 = SMat[:, 0, 0]
    squid.xaxis = squid.lin / flux0
    measdata.xaxis = np.linspace(-1+measdata.XPOS,
                                 1+measdata.XPOS,
                                 len(measdata.ydat))*measdata.XSC
    S11 = S11*10**(measdata.ATT / 20.0)
    S11 = addphase(S11, measdata.PHI)
    elem.S11 = S11
    return S11, measdata.ydat


def plotfig2(xaxis, xaxis2, S11, ydat):
    fig2 = plt.figure(2)
    fig2.clf()
    g1 = fig2.add_subplot(2, 1, 1)
    g1.hold(True)
    g1.plot(xaxis, abs(S11))
    g1.plot(xaxis2, abs(ydat))
    g1.hold(False)
    g2 = fig2.add_subplot(2, 1, 2)
    g2.hold(True)
    g2.plot(xaxis,  np.unwrap(np.angle(S11),  discont=pi))
    g2.plot(xaxis2, np.unwrap(np.angle(ydat), discont=pi))
    g2.hold(False)
    plt.draw()
    return


def plotfig5(xaxis, xaxis2, S11, ydat):
    fig5 = plt.figure(5)
    fig5.clf()
    g1 = fig5.add_subplot(2, 1, 1)
    g1.hold(True)
    g1.plot(xaxis, S11.real)
    g1.plot(xaxis2, ydat.real)
    g1.hold(False)
    g2 = fig5.add_subplot(2, 1, 2)
    g2.hold(True)
    g2.plot(xaxis,  S11.imag)
    g2.plot(xaxis2, ydat.imag)
    g2.hold(False)
    plt.draw()
    return


def find_nearest(someArray, value):
    '''
    Returns an index number (idx)
    at which the someArray.[idx] is closest to value
    '''
    idx = abs(someArray - value).argmin()
    return idx


def update(val):
    if elem.update is False:
        return
    f0 = sFreq.val * 1e9
    measdata.findex = find_nearest(measdata.freq, f0)
    measdata.ydat = measdata.D1complex[:, measdata.findex]
    squid.Ic = np.float64(sIc.val) * 1e-6
    squid.Cap = np.float64(sCap.val) * 1e-15
    squid.R = np.float64(sRsq.val) * 1e3
    elem.Z1 = np.float64(sZ1.val)
    elem.L1 = np.float64(sL1.val) * 1e-1
    elem.Z2 = np.float64(sZ2.val)
    elem.L2 = np.float64(sL2.val)
    elem.Z3 = np.float64(sZ3.val)
    elem.L3 = np.float64(sL3.val) * 1e-3
    elem.Z4 = np.float64(sZ4.val)
    measdata.XSC = np.float64(sXSC.val)
    measdata.XPOS = np.float64(sXPOS.val)
    measdata.ATT = np.float64(sATT.val)
    measdata.PHI = np.float64(sPHI.val)
    S11, ydat = getModelData(squid, elem, measdata)
    xaxis = squid.xaxis
    xaxis2 = measdata.xaxis
    plotfig2(xaxis, xaxis2, S11, ydat)
    plotfig5(xaxis, xaxis2, S11, ydat)
    return


def update2(val):
    Tstart = time()
    # Tell update not to run on each slider change
    elem.update = False
    sPHI.set_val(measdata.PHI)
    sATT.set_val(measdata.ATT)
    sIc.set_val(squid.Ic*1e6)
    sRsq.set_val(squid.R*1e-3)
    sCap.set_val(squid.Cap*1e15)
    sZ1.set_val(elem.Z1)
    sZ2.set_val(elem.Z2)
    elem.update = True
    sZ3.set_val(elem.Z3)
    # sL1.set_val(elem.L1*1e1)
    # sL2.set_val(elem.L2)
    # sL3.set_val(elem.L3*1e3)
    T = time()-Tstart
    print "sliders are updated ", T


def matchXaxis(val0):
    x2 = measdata.xaxis
    squid.start = x2[0]
    squid.stop = x2[-1]
    squid.update_lin(len(x2))
    squid.matchX = True
    update(0)
    return


def preFit(val0):
    S11, ydat = getModelData(squid, elem, measdata)
    zerofluxidx = find_nearest(measdata.xaxis, 0.0)
    t1 = np.unwrap(np.angle(ydat),  discont=pi)
    t2 = np.unwrap(np.angle(S11),  discont=pi)
    t3 = t1 - t2
    plt.figure(6)
    plt.clf()
    plt.plot(measdata.xaxis, (np.angle(S11) - np.angle(ydat)))
    plt.draw()
    measdata.PHI = measdata.PHI + t3[zerofluxidx]
    # t2 = 20*log10(np.abs(ydat)) - 20*log10(np.abs(S11))
    # measdata.ATT = measdata.ATT - t2[zerofluxidx]


def addphase(Data, Phi):
    Phi = Phi
    temp1 = np.empty(Data.shape, dtype='complex64')
    temp1 = (np.sin(Phi)*Data.real + np.cos(Phi)*Data.imag) * 1j
    temp1 = (np.cos(Phi)*Data.real - np.sin(Phi)*Data.imag) + temp1
    return temp1


def getfit():
    S11, ydat = getModelData(squid, elem, measdata)
    # S11 = S11*10**(measdata.ATT / 20.0)
    # S11 = addphase(S11, measdata.PHI)
    c = np.empty(len(S11)*2, dtype='float64')
    c[0::2] = S11.real
    c[1::2] = S11.imag
    return c


# def realimag(array):
#     return np.array([(x.real, x.imag) for x in array]).flatten()
# def residual(params, x, data=None):
#     ....
#     resid = calculate_complex_residual()
#     return realimag(resid)


def gta1(params, x, data):
    squid.R = params['R'].value
    squid.Cap = params['Cap'].value
    squid.Ic = params['Ic'].value
    elem.Z1 = params['Z1'].value
    # elem.Z3 = params['Z3'].value
    preFit(False)
    return getfit() - data


def fitcurve(val0):
    if squid.matchX is False:
        return
    ydat = measdata.D1complex[:, measdata.findex]
    data = np.empty(len(ydat)*2, dtype='float64')
    data[0::2] = ydat.real
    data[1::2] = ydat.imag
    preFit(False)
    xaxis3 = np.linspace(squid.start, squid.stop, (squid.pt*2))
    # Using standard curve_fit settings

    params = Parameters()
    params.add('Ic', value=squid.Ic, min=2.5e-6, max=4.5e-6)
    params.add('R', value=squid.R, min=1, max=1e5)
    params.add('Cap', value=squid.Cap, min=1e-15, max=1e-13)
    params.add('Z1', value=squid.Cap, min=25, max=100)
    params.add('Lbw', value=squid.Lbw, min=1e-90, max=1e-10)
    # params.add('Z3', value=squid.Cap, min=25, max=100)
    result = minimize(gta1, params, args=(xaxis3, data))
    print report_fit(result)

    squid.Ic = result.params['Ic'].value
    squid.R = result.params['R'].value
    squid.Cap = result.params['Cap'].value
    squid.Lbw = result.params['Lbw'].value
    elem.Z1 = result.params['Z1'].value
    # elem.Z3 = result.params['Z3'].value
    update2(0)
    preFit(True)

    # Calculate and plot residual there
    S11 = getfit()
    residual = data-S11
    plt.figure(4)
    plt.clf()
    plt.plot(xaxis3, residual)
    plt.draw()
    print abs(np.mean((residual*residual)))*1e8
    return


def xyFind(val):
    preFit(0)
    return


# --- Interface Buttons and Sliders ---
# --- Sliders Start ---
fig3 = plt.figure(3)
fig3.clear()
axcolor = 'lightgoldenrodyellow'
axATT = plt.axes([0.25, 0.80, 0.50, 0.03], axisbg=axcolor)
sATT = plt.Slider(axATT, 'Attenuation dBm', -90, -20.0,
                  valinit=-51.44, valfmt='%1.5f')
axPHI = plt.axes([0.25, 0.85, 0.50, 0.03], axisbg=axcolor)
sPHI = plt.Slider(axPHI, 'Phase offset', -np.pi, np.pi, valinit=0)
axXSC = plt.axes([0.25, 0.75, 0.50, 0.03], axisbg=axcolor)
sXSC = plt.Slider(axXSC, 'x-scale', 0.9, 1.1, valinit=1.04187, valfmt='%1.5f')
axXPOS = plt.axes([0.25, 0.70, 0.50, 0.03], axisbg=axcolor)
sXPOS = plt.Slider(axXPOS, 'x-pos', -0.5, 0.5, valinit=-0.49062, valfmt='%1.5f')
axZ4 = plt.axes([0.25, 0.60, 0.50, 0.03], axisbg=axcolor)
sZ4 = plt.Slider(axZ4, 'W.b. -> GND (Ohm)', 0.0001, 1.0, valinit=0.1)
axZ3 = plt.axes([0.25, 0.50, 0.50, 0.03], axisbg=axcolor)
sZ3 = plt.Slider(axZ3, 'Z3 (Ohm)', 0.0, 400.0, valinit=50)
axL3 = plt.axes([0.25, 0.55, 0.50, 0.03], axisbg=axcolor)
sL3 = plt.Slider(axL3, 'L3 (mm)', 0.0, 20.0, valinit=0.9)
axZ2 = plt.axes([0.25, 0.40, 0.50, 0.03], axisbg=axcolor)
sZ2 = plt.Slider(axZ2, 'Z2 (Ohm)', 0.0, 400.0, valinit=50)
axL2 = plt.axes([0.25, 0.45, 0.50, 0.03], axisbg=axcolor)
sL2 = plt.Slider(axL2, 'L2 (m)', 0.0, 1.0, valinit=0.3)
axZ1 = plt.axes([0.25, 0.30, 0.50, 0.03], axisbg=axcolor)
sZ1 = plt.Slider(axZ1, 'Z1 (Ohm)', 0.0, 900.0, valinit=50)
axL1 = plt.axes([0.25, 0.35, 0.50, 0.03], axisbg=axcolor)
sL1 = plt.Slider(axL1, 'L1 (m)', 0.01, 0.0167,
                 valinit=0.011625, valfmt='%1.10f')
axRsq = plt.axes([0.25, 0.25, 0.50, 0.03], axisbg=axcolor)
sRsq = plt.Slider(axRsq, 'Rsq (kOhm)', 0.01, 10.0, valinit=0.75)
axCap = plt.axes([0.25, 0.20, 0.50, 0.03], axisbg=axcolor)
sCap = plt.Slider(axCap, 'Cap (fF)', 0.01, 500.0, valinit=40)
axIc = plt.axes([0.25, 0.15, 0.50, 0.03], axisbg=axcolor)
sIc = plt.Slider(axIc, 'Ic (uA)', 0.1, 10.0, valinit=3.4)
axfreq = plt.axes([0.25, 0.1, 0.50, 0.03], axisbg=axcolor)
sFreq = plt.Slider(axfreq, 'Freq (GHz)', measdata.freq[0] / 1e9,
                   measdata.freq[-1] / 1e9, valinit=measdata.freq[0] / 1e9)

# --- Sliders Reactions ---
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

# --- Buttons
prFitxB = plt.axes([0.35, 0.025, 0.1, 0.04])
button5 = plt.Button(prFitxB, 'PreFit', color=axcolor, hovercolor='0.975')
button5.on_clicked(xyFind)

mXaxB = plt.axes([0.25, 0.025, 0.1, 0.04])
button4 = plt.Button(mXaxB, 'MatchX', color=axcolor, hovercolor='0.975')
button4.on_clicked(matchXaxis)

fitDaB = plt.axes([0.15, 0.025, 0.1, 0.04])
button3 = plt.Button(fitDaB, 'Fit', color=axcolor, hovercolor='0.975')
button3.on_clicked(fitcurve)

updatetax = plt.axes([0.05, 0.025, 0.1, 0.04])
button2 = plt.Button(updatetax, 'Update', color=axcolor, hovercolor='0.975')
button2.on_clicked(update2)

fig3.show()

# --- Interface End

update(0)
matchXaxis(0)
update(0)
