#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 16:55:40 2020

@author: val
"""
wave=6562.8
dwave1=1.2
dwave2=1.4

FWHM=46.7
dFWHM1=3.1
dFWHM2=2.6

###
dwave=(dwave1+dwave2)/2.0
dFWHM=(dFWHM1+dFWHM2)/2.0
Q=2*2.99792458e5*(FWHM/wave)
dQ=Q*((dwave/wave)**2 + (dFWHM/FWHM)**2)**0.5

print(Q)
print(dQ)


