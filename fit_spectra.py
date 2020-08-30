#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 23:16:15 2020

@author: val
"""
from astropy.io import fits
from astropy.convolution import Gaussian1DKernel, convolve
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal
import copy

import generic_peak

def renormalize(x,y,bounds):
    x_stripped = copy.deepcopy(x)
    y_stripped = copy.deepcopy(y)
    x_stripped=x_stripped.tolist()
    y_stripped=y_stripped.tolist()
    
    for i in range(len(x)):
        if x[i] > bounds[0] and x[i] < bounds[1]:
            x_stripped[i]=None
            y_stripped[i]=None
    x_stripped = [i for i in x_stripped if i]
    y_stripped = [i for i in y_stripped if i] 
    
    coeff,cov = np.polyfit(x_stripped, y_stripped, 1, cov=True)
    a,b=coeff
    err=np.sqrt(np.diag(cov))
    
    y_cor = y-(a*x+b)
    
    return y_cor,[a,b],err

fsdir="/Users/val/Desktop/"
fname="spec-0661-52163-0604"
fs=fits.open(fsdir+fname+".fits")

if fname[6] == "6": #0661
    z=0.198
elif fname[6] == "4": #0403
    z=0.312
else:
    z=0.0 
i=1

wavelength=fs[i].data.field('loglam')
wavelength=10**wavelength
wavelength=wavelength/(1+z)
flux=fs[i].data.field('flux')
err=(1/fs[i].data.field('ivar'))**0.5

gaussian = Gaussian1DKernel(stddev=7)
flux = convolve(flux,gaussian)
flux = scipy.signal.medfilt(flux,11)
    
plt.figure(figsize=(10,6))
plt.scatter(wavelength,flux,1.,color='k',marker='.')
#plt.plot(wavelength,flux,'k-')
#plt.errorbar(wavelength, flux, yerr=err, fmt=".k", capsize=0)
plt.xticks(np.arange(3000, 8500, step=500))
plt.xlabel("Angstroms")
plt.ylabel("1E-17 erg/cm^2/s/Ang")
plt.title(fname)
plt.show()

imin,imax=1175,2225

wavelength=wavelength[imin:imax]
flux=flux[imin:imax]
err=err[imin:imax]

#plt.figure(figsize=(10,6))
#plt.scatter(wavelength,flux,1.,color='k',marker='.')
#plt.plot(wavelength,flux,'k-')
#plt.errorbar(wavelength, flux, yerr=err, fmt=".k", capsize=0)
#plt.xticks(np.arange(6200, 6801, step=500))
#plt.xlabel("Angstroms")
#plt.ylabel("1E-17 erg/cm^2/s/Ang")
#plt.title(fname)
#plt.show()

x0,x1=4300,5100
fx_cor,p,line_err=renormalize(wavelength,flux,[x0,x1])

plt.figure(figsize=(10,6))
#plt.axvline(x=x0)
#plt.axvline(x=x1)
plt.scatter(wavelength,flux,1.,color='k',marker='.')
plt.plot(wavelength,p[0]*wavelength+p[1],'b')
plt.scatter(wavelength,fx_cor,1.,color='r',marker='.')
#plt.plot(wavelength,flux,'k-')
#plt.errorbar(wavelength, flux, yerr=err, fmt=".k", capsize=0)
plt.xticks(np.arange(4100, 5401, step=100))
plt.xlabel("Angstroms")
plt.ylabel("1E-17 erg/cm^2/s/Ang")
plt.title(fname)
plt.show()

################3
zoom=[475,850]

#x,y,y_err
x=wavelength[zoom[0]:zoom[1]]
y=fx_cor[zoom[0]:zoom[1]]
y_err=err[zoom[0]:zoom[1]]

#0661-Ha
#p1=[6565.0,15.0,20.0]
#p2=[6565.0,5.0,40.0]
#p3=[6725.0,2.0,10.0]

#0661-Hb
p1=[4860.0,3.0,10.0]
p2=[4850.0,2.7,39.0]
p3=[4895.0,0.9,10.0]

#p = parameters, p0 = initial guess
#p = f,a1,a2,3,u1,u2,u3,o1,o2,o3
p0=[0.1,p1[1],p2[1],p3[1],p1[0],p2[0],p3[0],p1[2],p2[2],p3[2]]

line = p[0]*x+p[1]
g1=generic_peak.gaussian(x,p0[1],p0[4],p0[7])
g2=generic_peak.gaussian(x,p0[2],p0[5],p0[8])
g3=generic_peak.gaussian(x,p0[3],p0[6],p0[9])
plt.figure(figsize=(10,6))
plt.plot(x, line+g1+g2+g3,color="b", alpha=0.8)
plt.plot(x, line+g1, color="r", alpha=0.8)
plt.plot(x, line+g2, color="r", alpha=0.8)
plt.plot(x, line+g3, color="r", alpha=0.8)
plt.scatter(wavelength[zoom[0]:zoom[1]],flux[zoom[0]:zoom[1]],1.,color='k',marker='.')
#plt.errorbar(wavelength[zoom[0]:zoom[1]], flux[zoom[0]:zoom[1]], yerr=err[zoom[0]:zoom[1]], fmt=".k",ms=2.0,elinewidth=1.0)
plt.legend(["initial","g1","g2","g3","data"])
plt.xlabel("Angstroms")
plt.ylabel("1E-17 erg/cm^2/s/Ang")
plt.title(fname)
plt.show()

#t = bounds on parameters
llim=[np.exp(-5.0),1e-8,1e-8,1e-8,4845,4835,4890,1e-8,1e-8,1e-8]
ulim=[np.exp(1.0),3.0*p0[1],3.0*p0[2],3.0*p0[3],4900,4853,4915,p0[7]+5.0,p0[8]+1.0,p0[9]+5.0]
t=[llim,ulim]

#emcee settings: w0 = # of walkers, w1 = # of iterations,
#w2 = parameter step size, #w3 = # of steps to discard
#w4 = # to flatten (thin)
w=[100,3000,1e-5,200,15]

#N = # of peaks to fit (1 to 3)
N=3

#plot = show total fit and data
#cor = show corner plots
#pos = show positions
#tau = calculate autocorrelation time

xfit,yfit,pfit,perr=generic_peak.fit_peaks(x,y,y_err,p0,t,w,N=N,plot=False,cor=False,pos=False,tau=False)

line = p[0]*xfit+p[1]
g1=line+generic_peak.gaussian(xfit,pfit[1],pfit[4],pfit[7])
g2=line+generic_peak.gaussian(xfit,pfit[2],pfit[5],pfit[8])
g3=line+generic_peak.gaussian(xfit,pfit[3],pfit[6],pfit[9])

x_pos,y_pos=x[0],17.0

textstr="continuum: \n"
textstr+='a = %.5f +/- %.5f \n' % (p[0],line_err[0])
textstr+='b = %.1f +/- %.1f \n' % (p[1],line_err[1])
textstr+="\ngaussian #1: \n"
textstr+='k = %.1f +/- (%.1f, %.1f) \n' % (pfit[1],perr[1][1],perr[0][1])
textstr+='u = %.1f +/- (%.1f, %.1f) \n' % (pfit[4],perr[1][4],perr[0][4])
textstr+='o = %.1f +/- (%.1f, %.1f) \n' % (pfit[7],perr[1][7],perr[0][7])
textstr+="\ngaussian #2: \n"
textstr+='k = %.1f +/- (%.1f, %.1f) \n' % (pfit[2],perr[1][2],perr[0][2])
textstr+='u = %.1f +/- (%.1f, %.1f) \n' % (pfit[5],perr[1][5],perr[0][5])
textstr+='o = %.1f +/- (%.1f, %.1f) \n' % (pfit[8],perr[1][8],perr[0][8])
textstr+="\ngaussian #3: \n"
textstr+='k = %.1f +/- (%.1f, %.1f) \n' % (pfit[3],perr[1][3],perr[0][3])
textstr+='u = %.1f +/- (%.1f, %.1f) \n' % (pfit[6],perr[1][6],perr[0][6])
textstr+='o = %.1f +/- (%.1f, %.1f)' % (pfit[9],perr[1][9],perr[0][9])


plt.figure(figsize=(20,12))
#plt.scatter(wavelength[zoom[0]:zoom[1]],flux[zoom[0]:zoom[1]],1.,color='k',marker='.')
plt.errorbar(wavelength[zoom[0]:zoom[1]], flux[zoom[0]:zoom[1]], yerr=err[zoom[0]:zoom[1]], fmt=".k",ms=3.0,elinewidth=1.0)
plt.plot(xfit, line+yfit,color="b", alpha=0.8)
plt.plot(xfit, g1, color="r", alpha=0.8)
plt.plot(xfit, g2, color="r", alpha=0.8)
plt.plot(xfit, g3, color="r", alpha=0.8)
plt.legend(["fit","g1","g2","g3","data"])
plt.text(x_pos, y_pos, textstr)
plt.xlabel("Angstroms")
plt.ylabel("1E-17 erg/cm^2/s/Ang")
plt.title(fname)
plt.savefig(fsdir+"0661_Hb_test.jpg",bbox_inches='tight')
plt.show()
