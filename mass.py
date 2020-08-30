#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 13:03:15 2019

@author: val
"""
import math
import numpy as np

def integral(k,u,o,kerr,uerr,oerr,a,b):
    t1=k*(np.sqrt(2.*np.pi)/2.)*o
    t1err=t1*np.sqrt((oerr/o)**2+(kerr/k)**2)
    
    x=(np.sqrt(2.)/2.)*(u-a)/o
    t2=math.erf(x)
    dx=x*np.sqrt((uerr/u)**2+(oerr/o)**2)
    t2err=(2/np.sqrt(np.pi))*(np.exp(-(x**2)))*dx
    
    x=(np.sqrt(2.)/2.)*(u-b)/o
    t3=math.erf(x)
    dx=x*np.sqrt((uerr/u)**2+(oerr/o)**2)
    t3err=(2/np.sqrt(np.pi))*(np.exp(-(x**2)))*dx
    
    I=t1*(t2-t3)
    Ierr=I*np.sqrt((t1err/t1)**2+((np.sqrt(t2err**2+t3err**2))/(t2-t3))**2)

    return I,Ierr #1E-17 erg/s/cm^2/Ang

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

def luminosity(z,k,u,o,kerr,uerr,oerr,a,b):
    DL=distance(z)
    I,Ierr = integral(k,u,o,kerr,uerr,oerr,a,b)
    L=4*np.pi*(DL**2)*I
    Lerr=4*np.pi*(DL**2)*Ierr
    return L,Lerr #1E-17 erg/s

def FWHM(u,o,uerr,oerr):
    F=2*2.99792458e5*(o/u)
    Ferr=F*np.sqrt((oerr/o)**2+(uerr/u)**2)
    return F,Ferr  #km/s

def log_mass_ratio(z,k,u,o,kerr,uerr,oerr,a,b,line):
    F,Ferr = FWHM(u,o,uerr,oerr)
    #print('FWHM:')
    #print(F)
    L,Lerr = luminosity(z,k,u,o,kerr,uerr,oerr,a,b)
    #print('Lum:')
    #print(L)
    
    if line == "ha":
        k1,k1err=6.3,0.0868
        k2,k2err=0.55,0.02
        k3,k3err=2.06,0.06
    elif line == "hb":
        k1,k1err=6.56,0.0241
        k2,k2err=0.56,0.02
        k3,k3err=2.0,0.0
    elif line == "mg2":
        k1,k1err=(1.70+3*2-2*0.63),0.07
        k2,k2err=0.63,0.00
        k3,k3err=2.0,0.0
    else:
        print("Emission line not supported.")
        print("Supported lines: ha, hb, mg2")
        return None
    
    R=k1+k2*np.log10(L*1e-59)+k3*np.log10(F*1e-3)
    
    Rerr = k1err**2
    Rerr+= k2**2 * (0.434*(Lerr/L))**2
    Rerr+= np.log10(L*1e-59)**2 * k2err**2
    Rerr+= k3**2 * (0.434*(Ferr/F))**2
    Rerr+= np.log10(F*1e-3)**2 * k3err**2
    Rerr = np.sqrt(Rerr)
    
    return R, Rerr #unitless 

# z=0.198    

# k,u,o=6.2,6562.8,46.7
# kerr,uerr,oerr=0.8,1.3,2.9
# a,b=6300.0,6800.0
# print(log_mass_ratio(z,k,u,o,kerr,uerr,oerr,a,b,"ha"))

# k,u,o=2.9,4851.2,39.4
# kerr,uerr,oerr=0.3,3.6,1.0
# a,b=4600.0,5100.0
# print(log_mass_ratio(z,k,u,o,kerr,uerr,oerr,a,b,"hb"))

# k,u,o=1.1,2799.1,38.1
# kerr,uerr,oerr=0.3,2.5,2.0
# a,b=2600,3000
# print(log_mass_ratio(z,k,u,o,kerr,uerr,oerr,a,b,"mg2"))

