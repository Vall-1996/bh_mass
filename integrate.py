#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 13:03:15 2020

@author: val
"""
import math
import numpy as np

def integral(k,u,o,a,b):
    t1=k*(np.sqrt(2.*np.pi)/2.)*o
    t2=math.erf((np.sqrt(2.)/2.)*(u-a)/o)
    t3=math.erf((np.sqrt(2.)/2.)*(u-b)/o)
    return t1*(t2-t3) #1E-17 erg/s/cm^2/Ang

def distance(z):
    H0,WM,WV=71.0,0.27,0.75
    
    h = H0/100.
    WR = 4.165E-5/(h*h)
    WK = 1-WM-WR-WV
    az = 1.0/(1+1.0*z)
    n=1000
    DTT = 0.0
    DCMR = 0.0
    
    for i in range(n):
        a = az+(1-az)*(i+0.5)/n
        adot = np.sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
        DTT = DTT + 1./adot
        DCMR = DCMR + 1./(a*adot)
    
    DCMR = (1.-az)*DCMR/n
    
    ratio = 1.00
    x = np.sqrt(np.abs(WK))*DCMR
    if x > 0.1:
        if WK > 0:
            ratio =  0.5*(np.exp(x)-np.exp(-x))/x 
        else:
            ratio = np.sin(x)/x
    else:
        y = x*x
        if WK < 0: y = -y
        ratio = 1. + y/6. + y*y/120.
    
    DCMT = ratio*DCMR
    DA = az*DCMT
    DL = DA/(az*az)
    DL_Mpc = (299792.458/H0)*DL
    
    return DL_Mpc*3.08568e24 #cm

def luminosity(z,k,u,o,a,b):
    DL=distance(z)
    I = integral(k,u,o,a,b)
    return 4*np.pi*(DL**2)*I #1E-17 erg/s

def FWHM(u,o):
    return 2*2.99792458e5*(o/u) #km/s

def log_mass_ratio(z,k,u,o,a,b,line):
    F = FWHM(u,o)
    L = luminosity(z,k,u,o,a,b)
    if line == "HA":
        k1=-32.33
        k2=0.55
        k3=2.06
    elif line == "HB":
        k1=-32.48
        k2=0.56
        k3=2.0
    else:
        print("Emission line not supported.")
        print("Supported lines: HA, HB")
        return None
    return k1+k2*np.log10(L)+k3*np.log10(F) #unitless 

z=0.198    
k,u,o=6.2,6562.8,46.7
a,b=6425.0,6725.0
print(log_mass_ratio(z,k,u,o,a,b,"HA"))





