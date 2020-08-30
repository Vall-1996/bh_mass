#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 20 17:20:18 2019

@author: val
"""
from astropy.convolution import Gaussian1DKernel, convolve
import scipy.signal
import matplotlib.pyplot as plt
import numpy as np

def renormalize(x,y,xb,plot=False):
    i1 = (np.abs(x - xb[0])).argmin()
    i2 = (np.abs(x - xb[1])).argmin()
    xs=np.concatenate((x[:i1],x[i2:]))
    ys=np.concatenate((y[:i1],y[i2:]))
    p,cov = np.polyfit(xs, ys, 1, cov=True)
    yc=y-linear(x,p[0],p[1])
    err=np.sqrt(np.diag(cov))
    
    if plot:
        plt.figure(figsize=(10,6))
        plt.axvline(x=xb[0],color="k",linestyle='-')
        plt.axvline(x=xb[1],color="k",linestyle='-')
        plt.scatter(x,y,1.,color='k',marker='.')
        plt.plot(x,linear(x,p[0],p[1]),'b')
        plt.scatter(x,yc,1.,color='r',marker='.')
        plt.xlabel("Angstroms")
        plt.ylabel("1E-17 erg/cm^2/s/Ang")
        plt.show()
        
    return yc,p,err

def trim(x,y,xb):
    i1 = (np.abs(x - xb[0])).argmin()
    i2 = (np.abs(x - xb[1])).argmin()
    return x[i1:i2],y[i1:i2]

def redshift(x,z):
    return x/(1+z)

def smooth(y,p=[7,11]):
    gaussian = Gaussian1DKernel(stddev=p[0])
    y = convolve(y,gaussian)
    return scipy.signal.medfilt(y,p[1])

def linear(x,a,b):
    return a*x+b

def gaussian(x,k,u,o):
    return k*np.exp(-0.5*((x-u)**2)/(o**2))

def display(x,y,err=[],p=[],p_err=[],gs=[],gs_nerr=[],gs_perr=[],t="",save="",big=False,cite=False,anchor=0):
    
    labels=["data"]
    
    if big:
        plt.figure(figsize=(20,12))
    else:
        plt.figure(figsize=(10,6))
        
    if len(err)>0:
        plt.errorbar(x,y,yerr=err,fmt=".k",ms=3.0,elinewidth=1.0,capsize=0)
    else:
        plt.scatter(x,y,1.,color='k',marker='.')
        
    if len(p)>0 and len(gs)<0:
        labels.append("line")
        plt.plot(x, linear(x,p[0],p[1]),color="b", alpha=0.8)
    if len(p)<0 and len(gs)>0:
        for i in range(len(gs)):
            labels.append("g"+str(i))
            plt.plot(x, gaussian(x,gs[i][0],gs[i][1],gs[i][2]), color="r", alpha=0.8)
    if len(p)>0 and len(gs)>0:
        for i in range(len(gs)):
            labels.append("g"+str(i))
            plt.plot(x, linear(x,p[0],p[1])+gaussian(x,gs[i][0],gs[i][1],gs[i][2]), color="r", alpha=0.8)
    if len(gs)>1:
        f=0
        labels.append("fit")
        for g in gs:
            f+=gaussian(x,g[0],g[1],g[2])
        if len(p)>0:
            f+=linear(x,p[0],p[1])
        plt.plot(x, f, color="b", alpha=0.8)
        
    if big and cite:
        textstr=""
        if len(p)>0:
            textstr="continuum: \n"
            textstr+='a = %.5f +/- %.5f \n' % (p[0],p_err[0])
            textstr+='b = %.1f +/- %.1f \n' % (p[1],p_err[1])
        if len(gs)>0:
            for i in range(len(gs)):
                textstr+="\ngaussian #"+str(i)+": \n"
                textstr+='k = %.1f +/- (%.1f, %.1f) \n' % (gs[i][0],gs_nerr[i][0],gs_perr[i][0])
                textstr+='u = %.1f +/- (%.1f, %.1f) \n' % (gs[i][1],gs_nerr[i][1],gs_perr[i][1])
                textstr+='o = %.1f +/- (%.1f, %.1f) \n' % (gs[i][2],gs_nerr[i][2],gs_perr[i][2])
        plt.text(x[0], anchor, textstr)
    
    if len(labels)>1:
        plt.legend(labels)
    plt.xlabel("Angstroms")
    plt.ylabel("1E-17 erg/cm^2/s/Ang")
    
    if len(t) > 0:
        plt.title(t)
        
    if len(save) > 0:
        import os
        fdir=os.getcwd()
        plt.savefig(fdir+"/"+save,bbox_inches='tight')
        
    plt.show()

