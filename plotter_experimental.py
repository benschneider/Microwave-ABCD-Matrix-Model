# -*- coding: utf-8 -*-
'''
Created on Fri Jul 10 15:02:34 2015
@author: benschneider
System:
SQUID at the end of a Transmission line.
ABCD-Matrix: M1, M2  each represent one element.
Ref: 'Microwave Engineering 3rd Edition' by David M. Pozar p. 185
'''
# from time import time
from numpy import pi, cos, abs, log10
import numpy as np
from parsers import dim, get_hdf5data
from ABCD import tline, sres, shunt, handler  # , terminator
# import matplotlib
# matplotlib.use('Qt4Agg')  # macosx, Qt4Agg, WX
# matplotlib.use('macosx')
# from scipy.optimize import curve_fit  # , leastsq
# from scipy.io import loadmat, savemat, whosmat #to save and load .mat (matlab)
from lmfit import minimize, Parameters, report_fit  # , Parameter
import matplotlib.pyplot as plt

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
elem.L2 = 0.04
elem.L3 = 0.01
elem.Z4 = 0.1           # Ohm; Wire bonds conductance to GND (-45dB isolation)
# m/s; approx. velocity in a coaxial 2/3 * speed of light
elem.v = 2.0e8
squid.Wb = 100e-12        # Aprox: Lw (nH) = 5.08x10-3 * L * (ln(4*L/D) - 1)
squid.LOOP = 1e-30
squid.ALP = 1.0
squid.f0 = 5e9

elem.updateOnChange = True     # Update when slider value is changed
elem.matchX = False


def get_sMatrix(elem, squid):
    '''
    This is where the assembly of the model happens using ABCD Matrixes
    SM represents the SQUID response
    M1 respresents 3 sections of Transmission lines (L1,Z1)(L2,Z2)(L3,Z3)
    Zsq is the impedance of the SQUID (flux) (vector)
    Assembly is done for each flux point and change in Zsq
    '''
    Zsq = get_Zsq(squid)
    b = 2.0 * pi * squid.f0 / 2.0e-8
    SM = np.zeros((len(Zsq), 2, 2)) * 1j  # complex matrix
    M1 = (tline(elem.Z1, b, elem.L1) *
          tline(elem.Z2, b, elem.L2) *
          tline(elem.Z3, b, elem.L3))  # transmission lines
    for ii, Zsq1 in enumerate(Zsq):
        M2 = sres(Zsq1) * shunt(elem.Z4)
        M4 = M1 * M2
        SM[ii] = elem.get_SM(M4)  # complex S-Matrix shape [2 2 fluxlength]
    return SM


def get_Zsq(squid):
    '''
    Calculates the Impedance of the SQUID.
    First calculates the inductance L  (flux)
    Then calulcates the Impedance (Cap, Freq, L)
    '''
    omega0 = 2.0 * pi * squid.f0
    flux = squid.lin
    # L = (flux0 / (2.0*pi*squid.Ic*abs(cos(pi*flux/squid.flux0))))
    # Ysq = (1.0/squid.R + 1j*omega0*squid.Cap - 1j/(omega0*(L + 1e-20)))
    # Alpha : Uneven junctions have an alpha other than 1.0
    # BC; alp = 0: Ic1=Ic, alp=2: Ic2 = Ic
    squid.Ic1 = squid.Ic*(1.0-squid.ALP/2.0)
    squid.Ic2 = squid.Ic - squid.Ic1
    print 'Ic1:', squid.Ic1, 'Ic2:', squid.Ic2
    L1 = (flux0 / (pi*squid.Ic1*abs(cos(pi*flux/squid.flux0))))
    L2 = (flux0 / (pi*squid.Ic2*abs(cos(pi*flux/squid.flux0))))
    Ysq1 = (0.5/squid.R + 0.5j*omega0*squid.Cap - 1j/(omega0*(L1 + 1e-90)))
    Ysq2 = (0.5/squid.R + 0.5j*omega0*squid.Cap - 1j/(omega0*(L2 + 1e-90)))
    Ysq = (Ysq1 + 1.0/(1.0/Ysq2 + 1j*(omega0*squid.LOOP)))
    return (1.0j*omega0*squid.Wb + 1.0 / Ysq)


def getModelData(elem, squid):
    '''
    Get S11 including Attenuation and Phase
    is using the full ABCD Model
    Phase and Attenuation is added (outside the ABCD model)
    '''
    squid.f0
    SMat = get_sMatrix(elem, squid)
    S11 = SMat[:, 0, 0]
    squid.xaxis = squid.lin / flux0
    measdata.xaxis = np.linspace(- 1 + measdata.XPOS, 1 + measdata.XPOS,
                                 len(measdata.ydat)) * measdata.XSC
    S11 = S11 * measdata.ATT  # * 10 ** (measdata.ATT / 20.0)
    S11 = addphase(S11, measdata.PHI)
    elem.S11 = S11


def plotfig2(xaxis, xaxis2, S11, ydat):
    fig2 = plt.figure(2)
    fig2.clf()
    g1 = fig2.add_subplot(2, 1, 1)
    g1.hold(True)
    g1.set_ylabel("Magnitude", fontsize=12)
    g1.plot(xaxis, abs(S11))
    g1.plot(xaxis2, abs(ydat))
    g1.axis('tight')
    g2 = fig2.add_subplot(2, 1, 2, sharex=g1)
    g2.hold(True)
    g2.set_ylabel("Phase rad.", fontsize=12)
    g2.set_xlabel("Flux/Flux0", fontsize=12)
    g2.plot(xaxis, np.unwrap(np.angle(S11), discont=pi))
    g2.plot(xaxis2, np.unwrap(np.angle(ydat), discont=pi))
    g2.axis('tight')
    plt.draw()
    return


def plotfig5(xaxis, xaxis2, S11, ydat):
    fig5 = plt.figure(5)
    fig5.clf()
    g1 = fig5.add_subplot(2, 1, 1)
    g1.hold(True)
    g1.set_ylabel("Real", fontsize=12)
    g1.plot(xaxis, S11.real)
    g1.plot(xaxis2, ydat.real)
    g1.axis('tight')
    g2 = fig5.add_subplot(2, 1, 2)
    g2.hold(True)
    g2.set_ylabel("Imag", fontsize=12)
    g2.set_xlabel("Flux/Flux0", fontsize=12)
    g2.plot(xaxis, S11.imag)
    g2.plot(xaxis2, ydat.imag)
    g2.axis('tight')
    plt.draw()
    return


def find_nearest(someArray, value):
    '''
    Returns an index number (idx)
    at which the someArray.[idx] is closest to value
    '''
    idx = abs(someArray - value).argmin()
    return idx


def update(val, doPlot=True):
    '''
    This Function updates the variables in elem and squid
    if the slider value is changed.
    doPlot can be set to false if no Plot is requested
    '''
    if elem.updateOnChange is False:
        # Do not run this function if no update is requested
        # This is done because all Sliders are set to autoupdate on change
        # Thus they will also call this when The sliders are updated.
        return
    squid.f0 = sFreq.val * 1e9
    measdata.findex = find_nearest(measdata.freq, squid.f0)
    measdata.ydat = measdata.D1complex[:, measdata.findex]
    squid.Ic = np.float64(sIc.val) * 1e-6
    squid.Cap = np.float64(sCap.val) * 1e-15
    squid.R = np.float64(sRsq.val) * 1e3
    squid.Wb = np.float64(sWb.val) * 1e-12
    elem.Z1 = np.float64(sZ1.val)
    elem.Z2 = np.float64(sZ2.val)
    elem.Z3 = np.float64(sZ3.val)
    elem.L1 = np.float64(sL1.val) * 1e-1
    elem.L2 = np.float64(sL2.val) * 1e-3
    elem.L3 = np.float64(sL3.val) * 1e-3
    elem.Z4 = np.float64(sZ4.val)
    measdata.XSC = np.float64(sXSC.val)
    measdata.XPOS = np.float64(sXPOS.val)
    measdata.ATT = np.float64(10 ** (sATT.val / 20))
    measdata.PHI = np.float64(sPHI.val)
    squid.ALP = np.float64(sALP.val)
    squid.LOOP = np.float64(sLOOP.val * 1e-9)
    getModelData(elem, squid)
    if doPlot is True:
        xaxis = squid.xaxis
        xaxis2 = measdata.xaxis
        plotfig2(xaxis, xaxis2, elem.S11, measdata.ydat)
        plotfig5(xaxis, xaxis2, elem.S11, measdata.ydat)
    return


def update2(val, update=True):
    '''
    It updates the sliders
    if Bool is True
    update() is also called
    if Bool is False sliders are only adjusted
    '''
    elem.updateOnChange = False  # Do not update figure on changes
    sPHI.set_val(measdata.PHI)
    sATT.set_val(20*log10(measdata.ATT))
    sIc.set_val(squid.Ic * 1e6)
    sRsq.set_val(squid.R * 1e-3)
    sCap.set_val(squid.Cap * 1e15)
    sWb.set_val(squid.Wb*1e12)
    sZ1.set_val(elem.Z1)
    sZ2.set_val(elem.Z2)
    sL2.set_val(elem.L2*1e3)
    sLOOP.set_val(squid.LOOP*1e9)
    sALP.set_val(squid.ALP)
    if update is True:
        elem.updateOnChange = True
    sZ3.set_val(elem.Z3)
    elem.updateOnChange = True


def matchXaxis(val0):
    x2 = measdata.xaxis
    squid.start = x2[0]
    squid.stop = x2[-1]
    squid.update_lin(len(x2))
    elem.matchX = True
    update(0)
    return


def preFit(val0):
    if elem.matchX is False:
        return
    getModelData(elem, squid)
    zerofluxidx = find_nearest(measdata.xaxis, 0.0)
    t2 = np.abs(measdata.ydat) / np.abs(elem.S11)
    t3 = (np.unwrap(np.angle(measdata.ydat), discont=pi) -
          np.unwrap(np.angle(elem.S11), discont=pi))
    measdata.ATT = measdata.ATT * t2[zerofluxidx]
    measdata.PHI = measdata.PHI + t3[zerofluxidx]


def addphase(Data, Phi):
    Phi = Phi
    temp1 = np.empty(Data.shape, dtype='complex64')
    temp1 = (np.sin(Phi) * Data.real + np.cos(Phi) * Data.imag) * 1j
    temp1 = (np.cos(Phi) * Data.real - np.sin(Phi) * Data.imag) + temp1
    return temp1


def getfit():
    getModelData(elem, squid)
    c = np.empty(len(elem.S11) * 2, dtype='float64')
    c[0::2] = elem.S11.real
    c[1::2] = elem.S11.imag
    return c


# def realimag(array):
#     return np.array([(x.real, x.imag) for x in array]).flatten()
# def residual(params, x, data=None):
#     ....
#     resid = calculate_complex_residual()
#     return realimag(resid)


def gta1(params, x, data):
    paramsToMem(params)
    print ('Ic:', squid.Ic, 'Wb:', squid.Wb, 'Cap:', squid.Cap)
    print ('Z1:', elem.Z1, 'Z2:', elem.Z2, 'Z3:', elem.Z3, 'L2:', elem.L2)
    preFit(False)
    return (getfit() - data)


def fitcurve(val0):
    if elem.matchX is False:
        return
    ydat = measdata.D1complex[:, measdata.findex]
    data = np.empty(len(ydat) * 2, dtype='float64')
    data[0::2] = ydat.real
    data[1::2] = ydat.imag
    preFit(False)
    xaxis3 = np.linspace(squid.start, squid.stop, (squid.pt * 2))
    # Using standard curve_fit settings

    # Define fitting parameters
    params = Parameters()
    params.add('CapfF', value=squid.Cap*1e15, vary=True, min=30, max=80)
    params.add('IcuA', value=squid.Ic*1e6, vary=True, min=3.0, max=4.0)
    params.add('WbpH', value=squid.Wb*1e12, vary=False, min=0, max=1000)
    params.add('LooppH', value=squid.LOOP*1e12, vary=True, min=0.0, max=100)
    params.add('alpha', value=squid.ALP, vary=True, min=0.98, max=1.02)
    params.add('R', value=squid.R, vary=True, min=1, max=20e3)
    params.add('Z1', value=elem.Z1, vary=True, min=40, max=60)
    params.add('Z2', value=elem.Z2, vary=True, min=40, max=60)
    params.add('Z3', value=elem.Z3, vary=True, min=40, max=60)
    params.add('L2', value=elem.L2, vary=False, min=0.00, max=0.09)

    # Do Fit
    result = minimize(gta1, params, args=(xaxis3, data))

    # Present results of fitting
    print report_fit(result)
    paramsToMem(result.params)
    update2(0)
    preFit(True)
    # Calculate and plot residual there
    S11 = getfit()
    residual = data - S11
    plt.figure(4)
    plt.clf()
    plt.plot(xaxis3, residual)
    plt.axis('tight')
    plt.draw()
    print 'Avg-sqr Residuals', abs(np.mean((residual * residual))) * 1e8
    return


def paramsToMem(params1):
    '''
    Help:
    How to do these things in a simple for loop ?
    something along the lines:
        for var in params
        varname is that of var
        obj.varname = var
    Maby this is easier to do in a Dict. format ?
    '''
    squid.Ic = params1['IcuA'].value*1e-6
    squid.Cap = params1['CapfF'].value*1e-15
    squid.Wb = params1['WbpH'].value*1e-12
    squid.LOOP = params1['LooppH'].value*1e-12
    squid.ALP = params1['alpha'].value
    squid.R = params1['R'].value
    elem.Z1 = params1['Z1'].value
    elem.Z2 = params1['Z2'].value
    elem.Z3 = params1['Z3'].value
    elem.L2 = params1['L2'].value


def preFitButton(val):
    preFit(0)
    update2(0)


def updateButton(val):
    if elem.updateOnChange is False:
        return
    update(val, doPlot=False)       # write sliders to mem
    preFit(0)                       # find new Phase & Att values
    update2(0, update=False)        # Update Sliders, elem.updateOnChange=False
    update(0)

# --- Interface Buttons and Sliders ---
# --- Sliders Start ---
'''
Some Interface stuff.
Also here it would be great to have some
automatic function which does this...
instead of writing each line,
something which detects what parameters are to be used and then
auto finds space and creates a slider for it.
'''
fig3 = plt.figure(3)
fig3.clear()
axcolor = 'lightgoldenrodyellow'
axATT = plt.axes([0.25, 0.90, 0.50, 0.02], axisbg=axcolor)
sATT = plt.Slider(axATT, 'Attenuation dBm', -90, -20.0,
                  valinit=-51.44, valfmt='%1.5f')
axPHI = plt.axes([0.25, 0.87, 0.50, 0.02], axisbg=axcolor)
sPHI = plt.Slider(axPHI, 'Phase offset', -np.pi, np.pi, valinit=0)
axXSC = plt.axes([0.25, 0.84, 0.50, 0.02], axisbg=axcolor)
sXSC = plt.Slider(axXSC, 'x-scale', 0.9, 1.1, valinit=1.04187, valfmt='%1.5f')
axXPOS = plt.axes([0.25, 0.81, 0.50, 0.02], axisbg=axcolor)
sXPOS = plt.Slider(axXPOS, 'x-pos', -0.5, 0.5, valinit=-0.49062, valfmt='%1.5f')

axWb = plt.axes([0.25, 0.49, 0.50, 0.02], axisbg=axcolor)
sWb = plt.Slider(axWb, 'WireBond Ind. pH', 0, 2000, valinit=0.0)
axALP = plt.axes([0.25, 0.46, 0.50, 0.02], axisbg=axcolor)
sALP = plt.Slider(axALP, 'Alpha', 0, 2, valinit=1)
axLOOP = plt.axes([0.25, 0.43, 0.50, 0.02], axisbg=axcolor)
sLOOP = plt.Slider(axLOOP, 'Loop Ind. nH', 0, 500, valinit=0.0)
axZ4 = plt.axes([0.25, 0.40, 0.50, 0.02], axisbg=axcolor)
sZ4 = plt.Slider(axZ4, 'Shrt-> GND (Ohm)', 0.0001, 1.0, valinit=0.1)
axZ3 = plt.axes([0.25, 0.37, 0.50, 0.02], axisbg=axcolor)
sZ3 = plt.Slider(axZ3, 'Z3 (Ohm)', 40.0, 60.0, valinit=50)
axL3 = plt.axes([0.25, 0.34, 0.50, 0.02], axisbg=axcolor)
sL3 = plt.Slider(axL3, 'L3 (mm)', 0.0, 20.0, valinit=0.0)
axZ2 = plt.axes([0.25, 0.31, 0.50, 0.02], axisbg=axcolor)
sZ2 = plt.Slider(axZ2, 'Z2 (Ohm)', 40.0, 60.0, valinit=50)
axL2 = plt.axes([0.25, 0.28, 0.50, 0.02], axisbg=axcolor)
sL2 = plt.Slider(axL2, 'L2 (mm)', 0.0, 80, valinit=0.0, valfmt='%1.5f')
axZ1 = plt.axes([0.25, 0.25, 0.50, 0.02], axisbg=axcolor)
sZ1 = plt.Slider(axZ1, 'Z1 (Ohm)', 40.0, 60.0, valinit=50)
axL1 = plt.axes([0.25, 0.22, 0.50, 0.02], axisbg=axcolor)
sL1 = plt.Slider(axL1, 'L1 (m)', 0.01, 0.0167, valinit=0.0, valfmt='%1.10f')
axRsq = plt.axes([0.25, 0.19, 0.50, 0.02], axisbg=axcolor)
sRsq = plt.Slider(axRsq, 'Rsq (kOhm)', 0.01, 10.0, valinit=0.75)
axCap = plt.axes([0.25, 0.16, 0.50, 0.02], axisbg=axcolor)
sCap = plt.Slider(axCap, 'Cap (fF)', 0, 400.0, valinit=60)
axIc = plt.axes([0.25, 0.13, 0.50, 0.02], axisbg=axcolor)
sIc = plt.Slider(axIc, 'Ic (uA)', 0.1, 10.0, valinit=3.4)
axfreq = plt.axes([0.25, 0.1, 0.50, 0.02], axisbg=axcolor)
sFreq = plt.Slider(axfreq, 'Freq (GHz)', measdata.freq[0] / 1e9,
                   measdata.freq[-1] / 1e9, valinit=measdata.freq[0] / 1e9)

# --- Sliders Reactions ---
sIc.on_changed(updateButton)
sCap.on_changed(updateButton)
sRsq.on_changed(updateButton)
sZ1.on_changed(updateButton)
sZ2.on_changed(updateButton)
sZ3.on_changed(updateButton)
sL1.on_changed(updateButton)
sL2.on_changed(updateButton)
sL3.on_changed(updateButton)
sZ4.on_changed(updateButton)
sXPOS.on_changed(updateButton)
sXSC.on_changed(updateButton)
sATT.on_changed(updateButton)
sPHI.on_changed(updateButton)
sFreq.on_changed(updateButton)
sLOOP.on_changed(updateButton)
sALP.on_changed(updateButton)
sWb.on_changed(updateButton)

# --- Buttons
prFitxB = plt.axes([0.35, 0.025, 0.1, 0.04])
button5 = plt.Button(prFitxB, 'PreFit', color=axcolor, hovercolor='0.975')
button5.on_clicked(preFitButton)

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

update(0, doPlot=False)
matchXaxis(0)
update(0, doPlot=False)
preFitButton(0)
