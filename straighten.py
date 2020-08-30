# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.

Created on Thu Mar 19 15:07:34 2020

@author: val
"""
import numpy as np
import matplotlib.pyplot as plt
import copy

#a,b = 0.3, 0.3
#x0,s,c = 50.,10.,2.0

def power(x,a,b):
    return a*(x**b)
def normal(x,x0,s,c):
    return c*np.exp(-0.5*((x-x0)**2)/s**2)

#x = np.linspace(10,100,1000)
#y = normal(x,x0,s,c) 
#y += power(x,a,b)

#plt.scatter(x,y)
#plt.show()

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
    
    a,b = np.polyfit(x_stripped, y_stripped, 1)
    y_cor = y-a*x+b
    m = min(y_cor)
    y_cor = y_cor-m
    
    return y_cor,[a,b,m]

#y_cor,p = renormalize(x,y,[20.0,80.0])

#plt.scatter(x,y_cor)
#plt.show()