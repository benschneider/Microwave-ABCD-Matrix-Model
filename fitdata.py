import numpy as np
from parsers import savemtx, loadmtx, make_header
from scipy.optimize import curve_fit  # , leastsq
# import scipy.optimize
from scipy.constants import Boltzmann as Kb
from scipy.constants import h, e  # , pi
from scipy.ndimage.filters import gaussian_filter1d
import matplotlib.pyplot as plt
from pymodelfit import FunctionModel1DAuto

# from pymodelfit import LinearModel


def fitcurve(measdata, squid, elem):
    '''
    measdata contains the data to be fittet for
    squid and elem are objects containing variables for the model.

    Fit Model > Data
    LSQ of (Data-Model)

    Model and Data contain complex numbers to be included in the fit.

    1. Ensure Model and Data -use Same number of points

    '''

    pass


class funcModel(FunctionModel1DAuto):
    def fitfun(self, f0, squid, elem):
        SMat = get_SMresponse(f0, squid, elem)
        S11 = SMat[:, 0, 0]
        return S11

