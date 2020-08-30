# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 12:16:58 2019
Depreciated tools for OzDES spectral fits files.
@author: div
"""

def integrate_voigt(pfit,lfit,wavelength_range,plot_values=False):
    from scipy.integrate import quad
    import numpy as np
    
    v1 = pfit[3]/np.sqrt(2*np.log(2))
    v2 = pfit[1]*(1-pfit[4])/(v1*np.sqrt(2*np.pi))
    v3 = 2*(v1**2)
    v4 = pfit[1]*pfit[4]*pfit[3]/np.pi
    voigt = lambda x: pfit[0]+(v2*np.exp(-((x-pfit[2])**2)/v3)+v4/((x-pfit[2])**2+pfit[3]**2))
    line = lambda x: lfit[0]*x+lfit[1]
    if plot_values:
        def voigt_t(x, b, A, u, o, a):
            g = o/np.sqrt(2*np.log(2))
            gauss = (1/(g*np.sqrt(2*np.pi)))*np.exp(-((x-u)**2)/(2*(g**2)))
            lorentz = (1/np.pi)*(o/((x-u)**2+o**2))
            return b+A*((1-a)*gauss+a*lorentz)
        xdata,ldata,pdata=[],[],[]
        for i in range(int(wavelength_range[0]-2),int(wavelength_range[1]+2)):
            xdata.append(i)
            pdata.append(voigt_t(i,pfit[0],pfit[1],pfit[2],pfit[3],pfit[4]))
            ldata.append(lfit[0]*i+lfit[1])
        return quad(voigt, wavelength_range[0], wavelength_range[1]), quad(line, wavelength_range[0], wavelength_range[1]),[xdata,pdata,ldata]
    return quad(voigt, wavelength_range[0], wavelength_range[1]), quad(line, wavelength_range[0], wavelength_range[1])

def integrate_simps(pfit,lfit,wavelength_range,step,plot=False):
    from scipy.integrate import simps
    import matplotlib.pyplot as plt
    import numpy as np
    def voigt(x, A, u, o, a):
        g = o/np.sqrt(2*np.log(2))
        gauss = (1/(g*np.sqrt(2*np.pi)))*np.exp(-((x-u)**2)/(2*(g**2)))
        lorentz = (1/np.pi)*(o/((x-u)**2+o**2))
        return A*((1-a)*gauss+a*lorentz)
    def voigts(x, pfit):
        sum_v=0.0
        if len(pfit)>2:
            for i in range(2,len(pfit),4):
                sum_v=sum_v+voigt(x,pfit[i],pfit[i+1],pfit[i+2],pfit[i+3])
        return pfit[0]*x+pfit[1]+sum_v
    def line(x, lfit):
        return lfit[0]*x+lfit[1]
    x = wavelength_range[0]
    xdata,pdata,ldata=[],[],[]
    while x<=wavelength_range[1]:
        xdata.append(x)
        pdata.append(voigts(x,pfit))
        ldata.append(line(x,lfit))
        x=x+step    
    if plot:
        plt.figure()
        plt.plot(xdata,pdata,'r-')
        plt.plot(xdata,ldata,'b-')
        plt.show()
    return (simps(pdata,xdata),0.0),(simps(ldata,xdata),0.0)

def find_intersects_voigt(lfit,pfit,wavelength_range,step):
    import numpy as np
    def voigt(x, b, A, u, o, a):
        g = o/np.sqrt(2*np.log(2))
        gauss = (1/(g*np.sqrt(2*np.pi)))*np.exp(-((x-u)**2)/(2*(g**2)))
        lorentz = (1/np.pi)*(o/((x-u)**2+o**2))
        voigt_x = b+A*((1-a)*gauss+a*lorentz)
        return voigt_x
    def line(x,a,b):
        return a*x+b
    x = wavelength_range[0]
    closeness=1000000
    intersect1 = -1
    while(x<=pfit[2]):
        peak = voigt(x,pfit[0],pfit[1],pfit[2],pfit[3],pfit[4])
        cont = line(x,lfit[0],lfit[1])
        if abs(peak-cont) < closeness:
            closeness = abs(peak-cont)
            intersect1=x
        x=x+step
    x=pfit[2]
    closeness=1000000
    intersect2 = -1
    while(x<=wavelength_range[1]):
        peak = voigt(x,pfit[0],pfit[1],pfit[2],pfit[3],pfit[4])
        cont = line(x,lfit[0],lfit[1])
        if abs(peak-cont) < closeness:
            closeness = abs(peak-cont)
            intersect2=x
        x=x+step
    return [intersect1,intersect2]

def fit_peaks(wavelength,flux,initial_fit,wavelength_range=[],slope=False,centers=False,plot=False,verbose=False):
    from scipy import optimize as opt
    import matplotlib.pyplot as plt
    import numpy as np
    def voigt(x, A, u, o, a):
        import numpy as np    
        g = o/np.sqrt(2*np.log(2))
        gauss = (1/(g*np.sqrt(2*np.pi)))*np.exp(-((x-u)**2)/(2*(g**2)))
        lorentz = (1/np.pi)*(o/((x-u)**2+o**2))
        voigt_x = A*((1-a)*gauss+a*lorentz)
        return voigt_x
    def function(x,*p):
        ret = p[0]*x+p[1]
        for i in range(2,len(p),4):
            ret = ret+voigt(x,p[i],p[i+1],p[i+2],p[i+3])
        return ret
    if len(wavelength_range) is 2:
        wavelength,flux=get_range(wavelength,flux,wavelength_range)
    initial_fit = [item for sublist in initial_fit for item in sublist]
    bounds_lower = (-0.00000001,-np.inf)
    bounds_upper = (0.00000001,np.inf)
    if slope:
        bounds_lower = (-np.inf,-np.inf)
        bounds_upper = (np.inf,np.inf)
    for i in range(2,len(initial_fit),4):
        bounds_lower = bounds_lower+(0.0,-np.inf,-np.inf,0.0)
        bounds_upper = bounds_upper+(np.inf,np.inf,np.inf,1.0)
    bounds = (bounds_lower,bounds_upper)
    popt, pcov = opt.curve_fit(function, wavelength, flux, initial_fit, bounds=bounds,maxfev=10000)
    perr = np.sqrt(np.diag(pcov))
    if verbose:
        print('m: '+str("{:.7f}".format(popt[0]))+' +/- '+str("{:.7f}".format(perr[0])))
        print('b: '+str("{:.7f}".format(popt[1]))+' +/- '+str("{:.7f}".format(perr[1])))
        print("")
        for i in range(2,len(initial_fit),4):
            print('A: '+str("{:.7f}".format(popt[i]))+' +/- '+str("{:.7f}".format(perr[i])))
            print('u: '+str("{:.7f}".format(popt[i+1]))+' +/- '+str("{:.7f}".format(perr[i+1])))
            print('o: '+str("{:.7f}".format(popt[i+2]))+' +/- '+str("{:.7f}".format(perr[i+2])))
            print('a: '+str("{:.7f}".format(popt[i+3]))+' +/- '+str("{:.7f}".format(perr[i+3])))
            print("")
    fdata = []
    for x in wavelength:
        fdata.append(function(x,*popt))
    if plot:
        plt.figure()
        plt.plot(wavelength, flux, 'k-')
        plt.plot(wavelength, fdata, 'r-')
        plt.show()
    if not centers:
        return wavelength,fdata,[popt,perr]
    flux_heights = []
    for i in range(2,len(popt),4):
        height_value = function(popt[i+1],*popt)
        width_value = 2*abs(popt[i+2])
        flux_heights.append([popt[i+1],height_value,width_value])
    return wavelength,fdata,[popt,perr],flux_heights