import matplotlib.pyplot as plt

''' This File contains the drawing of the boxes e.t.c'''

fig3 = plt.figure(3)
axcolor = 'lightgoldenrodyellow'
axfreq = plt.axes([0.25, 0.1, 0.50, 0.03], axisbg=axcolor)
axIc = plt.axes([0.25, 0.15, 0.50, 0.03], axisbg=axcolor)
axCap = plt.axes([0.25, 0.20, 0.50, 0.03], axisbg=axcolor)
axRsq = plt.axes([0.25, 0.25, 0.50, 0.03], axisbg=axcolor)
axZ1 = plt.axes([0.25, 0.30, 0.50, 0.03], axisbg=axcolor)
axL1 = plt.axes([0.25, 0.35, 0.50, 0.03], axisbg=axcolor)
axZ2 = plt.axes([0.25, 0.40, 0.50, 0.03], axisbg=axcolor)
axL2 = plt.axes([0.25, 0.45, 0.50, 0.03], axisbg=axcolor)
axZ3 = plt.axes([0.25, 0.50, 0.50, 0.03], axisbg=axcolor)
axL3 = plt.axes([0.25, 0.55, 0.50, 0.03], axisbg=axcolor)
axZ4 = plt.axes([0.25, 0.60, 0.50, 0.03], axisbg=axcolor)
axXPOS = plt.axes([0.25, 0.70, 0.50, 0.03], axisbg=axcolor)
axXSC = plt.axes([0.25, 0.75, 0.50, 0.03], axisbg=axcolor)
axATT = plt.axes([0.25, 0.80, 0.50, 0.03], axisbg=axcolor)
axPHI = plt.axes([0.25, 0.85, 0.50, 0.03], axisbg=axcolor)

sFreq = plt.Slider(axfreq,
                   'Freq (GHz)',
                   measdata.freq[0]/1e9,
                   measdata.freq[-1]/1e9,
                   valinit=measdata.freq[0]/1e9)
sIc = plt.Slider(axIc, 'Ic (uA)', 0.1, 10.0, valinit=3.4)
sCap = plt.Slider(axCap, 'Cap (fF)', 0.01, 500.0, valinit=40)
sRsq = plt.Slider(axRsq, 'Rsq (kOhm)', 0.01, 10.0, valinit=0.75)
sZ1 = plt.Slider(axZ1, 'Z1 (Ohm)', 0.0, 900.0, valinit=50)
sZ2 = plt.Slider(axZ2, 'Z2 (Ohm)', 0.0, 400.0, valinit=50)
sZ3 = plt.Slider(axZ3, 'Z3 (Ohm)', 0.0, 400.0, valinit=50)
sL1 = plt.Slider(axL1, 'L1 (m)', 0.0, 0.1, valinit=0.014062500000000006)
sL2 = plt.Slider(axL2, 'L2 (m)', 0.0, 1.0, valinit=0.3)
sL3 = plt.Slider(axL3, 'L3 (mm)', 0.0, 20.0, valinit=0.9)
sZ4 = plt.Slider(axZ4, 'W.b. -> GND (Ohm)', 0.0001, 1.0, valinit=0.1)
sXPOS = plt.Slider(axXPOS, 'x-pos', -1.0, 1.0, valinit=-0.49)
sXSC = plt.Slider(axXSC, 'x-scale', 0.5, 1.5, valinit=1.04)
sATT = plt.Slider(axATT, 'Attenuation dBm', -70, 0.0, valinit=-51.41)
sPHI = plt.Slider(axPHI, 'Phase offset', -360, 360.0, valinit=0)

fig3.show()


def find_nearest(someArray, value):
    idx = abs(someArray-value).argmin()
    return idx


def update(val):
    f0 = sFreq.val*1e9
    measdata.findex = find_nearest(measdata.freq, f0)
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
    SMat = get_SMresponse(f0, squid, elem)
    plotfig2(SMat, measdata)
    # plotfig4(SMat, measdata)

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

fitdata = plt.axes([0.3, 0.025, 0.1, 0.04])
button3 = plt.Button(fitdata, 'Fit', color=axcolor, hovercolor='0.975')
button3.on_clicked(fitcurve)

update(0)
