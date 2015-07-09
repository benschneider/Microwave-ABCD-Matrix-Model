# -*- coding: utf-8 -*-
'''
Created on Thu Jul  9 15:25:42 2015

@author: benschneider


Note:
Z = 1/Y
n = sqrt(er*ur) #refractive index n of the material
n ~ 3/2 #for a typical BNC cable
b = k/n = 2*pi/(wavelength*n)
c_light = 3e8 m/s
wavelength = c_light/frequency

2 Port circuits:
'''

import numpy as np

def impedance(Z):
    '''
    Impedance on one line
    ---|Z|---

    ---------
    M = [1 Z]
        [0 1]
    '''
    M = np.matrix([[1.0,Z],[0.0,1.0]])
    return M

def termination(Z):
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

def cable(Z,b,l):
    '''
    Cable with some length l, impedance Z, and adj. wave number b:
    ---------
     Z, b
    ---------
        l       #length

    M = [cos(bl)        i*Z*sin(b1)]
        [i*Y*sin(b1)    cos(bl)]
    '''
    Y = 1.0/Z
    M = np.matrix([[np.cos(b*l),1j*Z*np.sin(b*l)],[1j*Y*np.sin(b*l),np.cos(b*l)]])
    return M

def pulldown(Z):
    '''
    Pulldown or Pullup (1/Impedance to across both lines (i.e. to gnd))
    ----|----
        Y
    ----|----

    M = [1 0]
        [Y 1]
    '''
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