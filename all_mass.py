#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 22 16:02:15 2020

@author: val
"""
import numpy as np
import mass as mass


folder_path='/Users/val/Desktop/'
fileName = "fits.txt"
data = np.genfromtxt(folder_path+fileName, dtype=None)

ms,merrs=[],[]
for i in range(len(data)):
    line=data[i][0].decode('UTF-8')
    #print(line)
    z=data[i][1]
    k,u,o=data[i][2],data[i][4],data[i][6]
    kerr,uerr,oerr=data[i][3],data[i][5],data[i][7]
    a,b=u-200,u+200
    results = mass.log_mass_ratio(z,k,u,o,kerr,uerr,oerr,a,b,line)
    #print(results)
    #print("")
    ms.append(results[0])
    merrs.append(results[1])
final=[ms,merrs]
final=np.asarray(final)
final=np.transpose(final)

np.savetxt(folder_path+'results.txt', final , delimiter ='\t')  

fileName = "fits2.txt"
data = np.genfromtxt(folder_path+fileName, dtype=None)

ms,merrs=[],[]
for i in range(len(data)):
    line=data[i][0].decode('UTF-8')
    #print(line)
    z=data[i][1]
    k,u,o=data[i][2],data[i][4],data[i][6]
    kerr,uerr,oerr=data[i][3],data[i][5],data[i][7]
    a,b=u-200,u+200
    results = mass.log_mass_ratio(z,k,u,o,kerr,uerr,oerr,a,b,line)
    #print(results)
    #print("")
    ms.append(results[0])
    merrs.append(results[1])
final=[ms,merrs]
final=np.asarray(final)
final=np.transpose(final)

np.savetxt(folder_path+'results2.txt', final , delimiter ='\t')

