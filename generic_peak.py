#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 27 23:27:01 2019

@author: val
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as op
import emcee
import corner

def power_law(x,a,b,c):
    y=a*(x**b)+c
    return y

def gaussian(x,a,u,o):
    y=a*np.exp(-0.5*((x-u)**2)/(o**2))
    return y

#############################################################

def gen_data_3(N,f,p):
    u1,u2,u3,o1,o2,o3=p
    x = np.sort(10 * np.random.rand(N)+1)
    y_err = 0.05 + 0.05 * np.random.rand(N)
    y = gaussian(x,u1,o1)
    y += gaussian(x,u2,o2)
    y += gaussian(x,u3,o3)
    y += np.abs(f * y) * np.random.randn(N)
    y += y_err * np.random.randn(N)
    x_true = np.linspace(1, 11, 10*N)
    y_true = gaussian(x_true,u1,o1)
    y_true += gaussian(x_true,u2,o2)
    y_true += gaussian(x_true,u3,o3)
    return x,y,y_err,x_true,y_true

def gen_data_2(N,f,p):
    u1,u2,o1,o2=p
    x = np.sort(10 * np.random.rand(N)+1)
    y_err = 0.05 + 0.05 * np.random.rand(N)
    y = gaussian(x,u1,o1)
    y += gaussian(x,u2,o2)
    y += np.abs(f * y) * np.random.randn(N)
    y += y_err * np.random.randn(N)
    x_true = np.linspace(1, 11, 10*N)
    y_true = gaussian(x_true,u1,o1)
    y_true += gaussian(x_true,u2,o2)
    return x,y,y_err,x_true,y_true

def gen_data_1(N,f,p):
    u1,o1=p
    x = np.sort(10 * np.random.rand(N)+1)
    y_err = 0.05 + 0.05 * np.random.rand(N)
    y = gaussian(x,u1,o1)
    y += np.abs(f * y) * np.random.randn(N)
    y += y_err * np.random.randn(N)
    x_true = np.linspace(1, 11, 10*N)
    y_true = gaussian(x_true,u1,o1)
    return x,y,y_err,x_true,y_true

def gen_data(n,f,p,N=1,plot=True):
    
    if N == 1:
        x,y,y_err,x_true,y_true=gen_data_1(n,f,p)
    elif N == 2:
        x,y,y_err,x_true,y_true=gen_data_2(n,f,p)
    elif N == 3:
        x,y,y_err,x_true,y_true=gen_data_3(n,f,p) 
    else:
        print("Dimensionality supported: N=1,2,3.")
        return None, None, None, None, None
    
    if plot:
        plt.errorbar(x, y, yerr=y_err, fmt=".k", capsize=0)
        plt.plot(x_true, y_true, "r", alpha=0.8)
        plt.show()
        
    return x,y,y_err,x_true,y_true

#############################################################

def fit_peaks(x,y,yerr,p0,t,w,N=1,plot=True,cor=False,pos=False,tau=False):

    if N == 1:
        xfit,yfit,pfit,perr=single_fit(x,y,yerr,p0,t,w,plot,cor,pos,tau)
    elif N == 2:
        xfit,yfit,pfit,perr=double_fit(x,y,yerr,p0,t,w,plot,cor,pos,tau)
    elif N == 3:
        xfit,yfit,pfit,perr=triple_fit(x,y,yerr,p0,t,w,plot,cor,pos,tau)  
    else:
        print("Dimensionality supported: N=1,2,3.")
        return None, None, None, None
    
    return xfit,yfit,pfit,perr

def single_fit(x,y,y_err,p0,t,w,plot,cor,posi,tau):
    
    def func(x,a1,u1,o1):
        return gaussian(x,a1,u1,o1)
    
    def ln_like(p,x,y,y_err):
        lnf,a1,u1,o1 = p
        model = func(x,a1,u1,o1)
        noise_2 = y_err**2 + np.exp(2*lnf)*(model**2)
        items=(y-model)**2 / noise_2 + np.log(noise_2)
        return -0.5*np.sum(items)
    
    def ln_prior(p,t):
        lnf,a1,u1,o1 = p
        if np.log(t[0][0]) < lnf < np.log(t[1][0]) and t[0][1] < a1 < t[1][1] and t[0][2] < u1 < t[1][2] and t[0][3] < o1 < t[1][3]:
            return 0.0
        return -np.inf
    
    def ln_prob(p, t, x, y, y_err):
        lp = ln_prior(p,t)
        if not np.isfinite(lp):
            return -np.inf
        return lp + ln_like(p, x, y, y_err)
    
    f_0,a1_0,u1_0,o1_0=p0
    
    nll = lambda *args: -ln_like(*args)
    result = op.minimize(nll, [np.log(f_0),a1_0,u1_0,o1_0], args=(x, y, y_err))
    lnf_ml,a1_ml,u1_ml,o1_ml = result["x"]
    
    ndim, nwalkers = len(p0), w[0]
    pos = [result["x"] + w[2]*np.random.randn(ndim) for i in range(nwalkers)]
    
    sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_prob, args=(t, x, y, y_err))
    sampler.run_mcmc(pos, w[1], progress=True)
    samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
    
    if cor:
        fig = corner.corner(samples,labels=["$\ln\,f$", "$a1$", "$u1$", "$o1$"])
        plt.show(fig)
        
    if posi:    
        fig, axes = plt.subplots(len(p0), figsize=(10, 6), sharex=True)
        samples = sampler.get_chain()
        labels = ["lf ", "a1", "u1", "o1"]
        for i in range(ndim):
            ax = axes[i]
            ax.plot(samples[:, :, i], "k", alpha=0.3)
            ax.set_xlim(0, len(samples))
            ax.set_ylabel(labels[i])
            ax.yaxis.set_label_coords(-0.1, 0.5)
        axes[-1].set_xlabel("step number")
        plt.show()
        
    if tau:
        dt = np.mean(sampler.get_autocorr_time())
        print("AC time = "+'{:.3f}'.format(dt))
        w[3]=int(2.5*dt)
        w[4]=int(0.375*dt)
    
    xfit=np.linspace(x[0], x[-1], 100)
    pfit,perr=[],[[],[]]
    
    flat_samples = sampler.get_chain(discard=w[3], thin=15, flat=True)
    labels = ["f ", "a1", "u1", "o1"]
    
    logf = np.percentile(flat_samples[:, 0], [16, 50, 84])
    logq = np.diff(logf)
    pfit.append(np.exp(logf[1]))
    perr[0].append(np.exp(logf[1])*logq[0])
    perr[1].append(np.exp(logf[1])*logq[1])
    
    print(labels[0]+" = "+'{:.3f}'.format(np.exp(logf[1]))+" - "+'{:.3f}'.format(np.exp(logf[1])*logq[0])+" + "+'{:.3f}'.format(np.exp(logf[1])*logq[1]))

    for i in range(1,ndim):
        mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        pfit.append(mcmc[1])
        perr[0].append(q[0])
        perr[1].append(q[1])
        print(labels[i]+" = "+'{:.3f}'.format(mcmc[1])+" - "+'{:.3f}'.format(q[0])+" + "+'{:.3f}'.format(q[1]))

    yfit=func(xfit,pfit[1],pfit[2],pfit[3])
    y0=func(xfit,a1_0,u1_0,o1_0)
    
    if plot:
        plt.figure(figsize=(10,6))
        plt.plot(xfit, y0, color="b", alpha=0.8)
        plt.plot(xfit, yfit, color="r", alpha=0.8)
        plt.errorbar(x, y, yerr=y_err, fmt=".k",ms=1.0,elinewidth=1.0)
        plt.legend(["guess","fit","data"])
        plt.show()
    
    return xfit,yfit,pfit,perr

def double_fit(x,y,y_err,p0,t,w,plot,cor,posi,tau):
    
    def func(x,a1,a2,u1,u2,o1,o2):
        y = gaussian(x,a1,u1,o1)
        y += gaussian(x,a2,u2,o2)
        return y
    
    def ln_like(p,x,y,y_err):
        lnf,a1,a2,u1,u2,o1,o2 = p
        model = func(x,a1,a2,u1,u2,o1,o2)
        noise_2 = y_err**2 + np.exp(2*lnf)*(model**2)
        items=(y-model)**2 / noise_2 + np.log(noise_2)
        return -0.5*np.sum(items)
    
    def ln_prior(p,t):
        lnf,a1,a2,u1,u2,o1,o2 = p
        if np.log(t[0][0]) < lnf < np.log(t[1][0]) and t[0][1] < a1 < t[1][1] and t[0][2] < a2 < t[1][2] and t[0][3] < u1 < t[1][3] and t[0][4] < u2 < t[1][4] and t[0][5] < o1 < t[1][5] and t[0][6] < o2 < t[1][6]:
            return 0.0
        return -np.inf
    
    def ln_prob(p, t, x, y, y_err):
        lp = ln_prior(p,t)
        if not np.isfinite(lp):
            return -np.inf
        return lp + ln_like(p, x, y, y_err)
    
    f_0,a1_0,a2_0,u1_0,u2_0,o1_0,o2_0=p0
    
    nll = lambda *args: -ln_like(*args)
    result = op.minimize(nll, [np.log(f_0),a1_0,a2_0,u1_0,u2_0,o1_0,o2_0], args=(x, y, y_err))
    lnf_ml,a1_ml,a2_ml,u1_ml,u2_ml,o1_ml,o2_ml = result["x"]
    
    ndim, nwalkers = len(p0), w[0]
    pos = [result["x"] + w[2]*np.random.randn(ndim) for i in range(nwalkers)]
    
    sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_prob, args=(t, x, y, y_err))
    sampler.run_mcmc(pos, w[1], progress=True)
    samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
    
    if cor:
        fig = corner.corner(samples,labels=["$\ln\,f$", "a1", "a2","$u1$", "$u2$", "$o1$", "$o2$"])
        plt.show(fig)
        
    if posi:    
        fig, axes = plt.subplots(len(p0), figsize=(10, 8), sharex=True)
        samples = sampler.get_chain()
        labels = ["lf ", "a1", "a2", "u1", "u2", "o1", "o2"]
        for i in range(ndim):
            ax = axes[i]
            ax.plot(samples[:, :, i], "k", alpha=0.3)
            ax.set_xlim(0, len(samples))
            ax.set_ylabel(labels[i])
            ax.yaxis.set_label_coords(-0.1, 0.5)
        axes[-1].set_xlabel("step number")
        plt.show()
        
    if tau:
        dt = np.mean(sampler.get_autocorr_time())
        print("AC time = "+'{:.3f}'.format(dt))
        w[3]=int(2.5*dt)
        w[4]=int(0.375*dt)
    
    xfit=np.linspace(x[0], x[-1], 100)
    pfit,perr=[],[[],[]]
    
    flat_samples = sampler.get_chain(discard=w[3], thin=w[4], flat=True)
    labels = ["f ", "a1", "a2", "u1", "u2", "o1", "o2"]
    
    logf = np.percentile(flat_samples[:, 0], [16, 50, 84])
    logq = np.diff(logf)
    pfit.append(np.exp(logf[1]))
    perr[0].append(np.exp(logf[1])*logq[0])
    perr[1].append(np.exp(logf[1])*logq[1])
    
    print(labels[0]+" = "+'{:.3f}'.format(np.exp(logf[1]))+" - "+'{:.3f}'.format(np.exp(logf[1])*logq[0])+" + "+'{:.3f}'.format(np.exp(logf[1])*logq[1]))

    for i in range(1,ndim):
        mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        pfit.append(mcmc[1])
        perr[0].append(q[0])
        perr[1].append(q[1])
        print(labels[i]+" = "+'{:.3f}'.format(mcmc[1])+" - "+'{:.3f}'.format(q[0])+" + "+'{:.3f}'.format(q[1]))

    yfit=func(xfit,pfit[1],pfit[2],pfit[3],pfit[4],pfit[5],pfit[6])
    y0=func(xfit,a1_0,a2_0,u1_0,u2_0,o1_0,o2_0)
    
    if plot:
        plt.figure(figsize=(10,6))
        plt.plot(xfit, y0, color="b", alpha=0.8)
        plt.plot(xfit, yfit, color="r", alpha=0.8)
        plt.errorbar(x, y, yerr=y_err, fmt=".k",ms=1.0,elinewidth=1.0)
        plt.legend(["guess","fit","data"])
        plt.show()
    
    return xfit,yfit,pfit,perr

def triple_fit(x,y,y_err,p0,t,w,plot,cor,posi,tau):
    
    def func(x,a1,a2,a3,u1,u2,u3,o1,o2,o3):
        y = gaussian(x,a1,u1,o1)
        y += gaussian(x,a2,u2,o2)
        y += gaussian(x,a3,u3,o3)
        return y
    
    def ln_like(p,x,y,y_err):
        lnf,a1,a2,a3,u1,u2,u3,o1,o2,o3 = p
        model = func(x,a1,a2,a3,u1,u2,u3,o1,o2,o3)
        noise_2 = y_err**2 + np.exp(2*lnf)*(model**2)
        items=(y-model)**2 / noise_2 + np.log(noise_2)
        return -0.5*np.sum(items)
    
    def ln_prior(p,t):
        lnf,a1,a2,a3,u1,u2,u3,o1,o2,o3 = p
        if np.log(t[0][0]) < lnf < np.log(t[1][0]) and t[0][1] < a1 < t[1][1] and t[0][2] < a2 < t[1][2] and t[0][3] < a3 < t[1][3] and t[0][4] < u1 < t[1][4] and t[0][5] < u2 < t[1][5] and t[0][6] < u3 < t[1][6] and t[0][7] < o1 < t[1][7] and t[0][8] < o2 < t[1][8] and t[0][9] < o3 < t[1][9]:
            return 0.0
        return -np.inf
    
    def ln_prob(p, t, x, y, y_err):
        lp = ln_prior(p,t)
        if not np.isfinite(lp):
            return -np.inf
        return lp + ln_like(p, x, y, y_err)
    
    f_0,a1_0,a2_0,a3_0,u1_0,u2_0,u3_0,o1_0,o2_0,o3_0=p0
    
    nll = lambda *args: -ln_like(*args)
    result = op.minimize(nll, [np.log(f_0),a1_0,a2_0,a3_0,u1_0,u2_0,u3_0,o1_0,o2_0,o3_0], args=(x, y, y_err))
    lnf_ml,a1_ml,a2_ml,a3_ml,u1_ml,u2_ml,u3_ml,o1_ml,o2_ml,o3_ml = result["x"]
    
    ndim, nwalkers = len(p0), w[0]
    pos = [result["x"] + w[2]*np.random.randn(ndim) for i in range(nwalkers)]
    
    sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_prob, args=(t, x, y, y_err))
    sampler.run_mcmc(pos, w[1], progress=True)
    samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
    
    if cor:
        fig = corner.corner(samples,labels=["$\ln\,f$", "$a1$", "$a2$", "$a3$", "$u1$", "$u2$","$u3$", "$o1$", "$o2$", "$o3$"])
        plt.show(fig)
    
    if posi:    
        fig, axes = plt.subplots(len(p0), figsize=(10, 10), sharex=True)
        samples = sampler.get_chain()
        labels = ["lf ", "a1", "a2", "a3", "u1", "u2", "u3", "o1", "o2", "o3"]
        for i in range(ndim):
            ax = axes[i]
            ax.plot(samples[:, :, i], "k", alpha=0.3)
            ax.set_xlim(0, len(samples))
            ax.set_ylabel(labels[i])
            ax.yaxis.set_label_coords(-0.1, 0.5)
        axes[-1].set_xlabel("step number")
        plt.show()
        
    if tau:
        dt = np.mean(sampler.get_autocorr_time())
        print("AC time = "+'{:.3f}'.format(dt))
        w[3]=int(2.5*dt)
        w[4]=int(0.375*dt)
    
    xfit=np.linspace(x[0], x[-1], 1000)
    pfit,perr=[],[[],[]]
    
    flat_samples = sampler.get_chain(discard=w[3], thin=w[4], flat=True)
    labels = ["f ", "a1", "a2", "a3", "u1", "u2", "u3", "o1", "o2", "o3"]
    
    logf = np.percentile(flat_samples[:, 0], [16, 50, 84])
    logq = np.diff(logf)
    pfit.append(np.exp(logf[1]))
    perr[0].append(np.exp(logf[1])*logq[0])
    perr[1].append(np.exp(logf[1])*logq[1])
    
    print(labels[0]+" = "+'{:.3f}'.format(np.exp(logf[1]))+" - "+'{:.3f}'.format(np.exp(logf[1])*logq[0])+" + "+'{:.3f}'.format(np.exp(logf[1])*logq[1]))

    for i in range(1,ndim):
        mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
        q = np.diff(mcmc)
        pfit.append(mcmc[1])
        perr[0].append(q[0])
        perr[1].append(q[1])
        print(labels[i]+" = "+'{:.3f}'.format(mcmc[1])+" - "+'{:.3f}'.format(q[0])+" + "+'{:.3f}'.format(q[1]))

    yfit=func(xfit,pfit[1],pfit[2],pfit[3],pfit[4],pfit[5],pfit[6],pfit[7],pfit[8],pfit[9])
    y0=func(xfit,a1_0,a2_0,a3_0,u1_0,u2_0,u3_0,o1_0,o2_0,o3_0)
    
    if plot:
        plt.figure(figsize=(10,6))
        plt.plot(xfit, y0, color="b", alpha=0.8)
        plt.plot(xfit, yfit, color="r", alpha=0.8)
        plt.errorbar(x, y, yerr=y_err, fmt=".k",ms=1.0,elinewidth=1.0)
        plt.legend(["guess","fit","data"])
        plt.show()
        
    return xfit,yfit,pfit,perr

#############################################################

#u1_true,o1_true = 6.0,0.25
#u2_true,o2_true = 4.0,0.5
#u3_true,o3_true = 5.0,0.35
#f_true = 0.1

#N = 50
#p3=[u1_true,u2_true,u3_true,o1_true,o2_true,o3_true]
#x,y,y_err,x_true,y_true = gen_data(N,f_true,p3,3,plot=False)

#############################################################

# #x,y,y_err

# #p = parameters, p0 = initial guess
# #p = f,u1,u2,u3,o1,o2,o3
# p0=[0.15,5.9,4.1,5.2,0.3,0.6,0.6]

# #t = bounds on parameters
# llim=[np.exp(-5.0),x[0],x[0],x[0],1e-8,1e-8,1e-8]
# ulim=[np.exp(1.0),x[-1],x[-1],x[-1],1,1,1]
# t=[llim,ulim]

# #emcee settings: w0 = # of walkers, w1 = # of iterations,
# #w2 = parameter step size, #w3 = # of steps to discard
# #w4 = # to flatten (thin)
# w=[250,6000,1e-5,200,15]

# #N = # of peaks to fit (1 to 3)
# N=3

# #plot = show total fit and data
# #cor = show corner plots
# #pos = show positions
# #tau = calculate autocorrelation time

# xfit,yfit,pfit,perr=fit_peaks(x,y,y_err,p0,t,w,N=N,plot=True,cor=False,pos=True,tau=True)


