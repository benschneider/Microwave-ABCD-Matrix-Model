# -*- coding: utf-8 -*-
'''
Created on Thu Jul  9 15:25:42 2015
@author: benschneider
'''

import numpy as np

#class components():
'''
Z = 1/Y
n = sqrt(er*ur) #refractive index n of the material
n ~ 3/2 #for a typical BNC cable
b = k/n = 2*pi/(wavelength*n)
c_light = 3e8 m/s
wavelength = c_light/frequency

2 Port circuits:
Ref: 'Microwave Engineering 3rd Edition' by David M. Pozar p. 185
'''

def __init__(self):
    pass

#__all__ = ["sres", "terminator", "tline", "shunt", "coupling", "handler"]

def sres(Z):
    '''
    Series resistor on one line
    ---|Z|---

    ---------
    M = [1 Z]
        [0 1]
    '''
    M = np.matrix([[1.0,Z],[0.0,1.0]])
    return M

def terminator(Z):
    '''
    Termination (Impedance on line & Port 2 shorted)
    --|Z|--|
           |
    -------|
    M = [0 Z]
        [0 1]
    Det = 0 ! #non reciprocal network!
    '''
    M = np.matrix([[0.0,Z],[0.0,1.0]])
    return M

def tline(Z,b,l):
    '''
    Transmission line with: length l, impedance Z, and adj. wave number b:
    ---------
     Z, b
    ---------
        l       #length

    M = [cos(bl)        i*Z*sin(b1)]
        [i*Y*sin(b1)    cos(bl)]
    '''
    if Z == 0:
        Z = 1e-100
    Y = 1.0/Z
    M = np.matrix([[np.cos(b*l),1j*Z*np.sin(b*l)],[1j*Y*np.sin(b*l),np.cos(b*l)]])
    return M

def shunt(Z):
    '''
    Shunt resistor (1/Impedance to across both lines (i.e. to gnd))
    ----|----
        Y
    ----|----

    M = [1 0]
        [Y 1]
    '''
    if Z == 0:
        Z = 1e-100

    Y = 1.0/Z
    M = np.matrix([[1.0,0.0],[Y,1.0]])
    return M

def coupling(N):
    '''
    # Coupling
       N:1
    ---} {-----
       { }
    ---} {-----

    M = [N   0]
        [0   1/N]
    '''
    M = np.matrix([[N,0.0],[0,1.0/N]])
    return M


class handler():
    '''
    Handling data recording
    also used to convert ABCD element results
    '''
    def __init__(self, name = 'mag/phase' ,start = 0, stop = 7, pt = 8):
        self.name = name     #ufo: unknown fried object
        self.start = start
        self.stop = stop
        self.pt = int(pt)
        self.Z0 = 50
        self.lin = np.linspace(self.start,self.stop,self.pt)

    def update_lin(self,pt):
        self.pt = int(pt)
        self.lin = np.linspace(self.start,self.stop,self.pt)

    def prepare_data_save(self, dim_1, dim_2, dim_3):
        self._SMat = np.zeros((dim_3.pt, dim_2.pt, dim_1.pt))

    def record_SM(self,ABCD, jj,ii):
        '''
        Stores the results of the S Matrix into '_SMat'
        '''
        S = self.get_SM(ABCD)
        self._SMat[0,jj,ii]     = np.abs(S[0,0])
        self._SMat[1,jj,ii]     = np.angle(S[0,0])
        self._SMat[2,jj,ii]     = np.abs(S[0,1])
        self._SMat[3,jj,ii]     = np.angle(S[0,1])
        self._SMat[4,jj,ii]     = np.abs(S[1,0])
        self._SMat[5,jj,ii]     = np.angle(S[1,0])
        self._SMat[6,jj,ii]     = np.abs(S[1,1])
        self._SMat[7,jj,ii]     = np.angle(S[1,1])

    def record_ZL(self,Z,L, jj,ii):
        '''
        Stores the results of the S Matrix into '_SMat'
        '''
        self._SMat[8,jj,ii] = np.abs(Z)
        self._SMat[9,jj,ii] = np.angle(Z)
        self._SMat[10,jj,ii] = np.abs(L)

    def unwrap_SM(self,jj):
        #a better unwrap function is still in thinking...
        self._SMat[1,jj]     = np.unwrap(self._SMat[1,jj])
        self._SMat[3,jj]     = np.unwrap(self._SMat[3,jj])
        self._SMat[5,jj]     = np.unwrap(self._SMat[5,jj])
        self._SMat[7,jj]     = np.unwrap(self._SMat[6,jj])

    def get_SM(self, M):
        Z0 = self.Z0
        '''
        Returns S Matrix
        [S11,S12]
        [S21,S22]
        '''
        A = M[0,0]
        B = M[0,1]
        C = M[1,0]
        D = M[1,1]
        S11 = ( A+B/Z0-C*Z0-D ) / ( A+B/Z0+C*Z0+D )
        S12 = 2.0*( A*D-B*C ) / ( A+B/Z0+C*Z0+D )
        S21 = 2.0 / ( A+B/Z0+C*Z0+D )
        S22 = ( -A+B/Z0-C*Z0+D ) / ( A+B/Z0+C*Z0+D )
        return np.matrix([[S11,S12],[S21,S22]])