#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 16:02:15 2019

@author: val
"""

#from astropy.io import fits
#from astropy.table import Table, Column
#from astropy.coordinates import SkyCoord, Angle
#import astropy.units as u
#from astroquery.irsa import Irsa
#from openpyxl import Workbook 
#from openpyxl import load_workbook
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.patches as mpatches
import ozdes_tools as oz

folder_path='/Users/val/Desktop/university/Research/SpARCS_data/'
#filename='SVA1_COADD-2971189603.fits'
#filename='SVA1_COADD-2938301050.fits'
#filename='SVA1_COADD-2938396262.fits'
#filename='SVA1_COADD-2940250491.fits'

##### Determining redshift of a single spectra.

files = oz.summon_files(folder_path)

f=2
m=[0]
i=0
wavelength,flux,z = oz.get_spectra(folder_path, files[f], i=i, m=m,lines=True,plot=True)
#print(files[f])



# #initial_fit=[70,2800,50,1.5]
# #intersect_limit = [2750,2850]
# #c,d=[1.70,0.07],[0.63,0.00]
# initial_fit=[70,4863,50,1.5]
# intersect_limit = [4760,4960]
# c,d=[1.63,0.04],[0.49,0.03]
# #initial_fit=[70,6565,50,1.5]
# #intersect_limit = [6465,6665]
# #c,d=[1.0,1.0],[1.0,1.0]
# ks = [15,10,2,2]
# plots=[False,True,False]
# zoom = False
# DL=-1
# oz.autorun_fit(folder_path,files[f],initial_fit,m,ks,intersect_limit,plots,zoom,DL,c,d)



# #for filename in files:
# #    wavelength,flux,z=oz.get_spectra(folder_path, filename, i=0, m=[], mask_break=True,lines=True,save=True)
#     #xs,dxs,dzs=[2799.117,4862.68,6564.61],[40,70,100],[0.01,0.02,0.03]
#     #oz.plot_emissions(folder_path,filename,wavelength,flux,z,xs=xs,dxs=dxs,dzs=dzs)


# #mb,p,ls=True,False,True
# #wavelength,flux,z=oz.get_spectra(folder_path, filename, i=0, m=[], mask_break=mb, plot=p, lines=ls)
# #wavelength_1,flux_1=oz.adjust_spectra(wavelength,flux,z,dz=0.05)
# #wavelength_2,flux_2=oz.adjust_spectra(wavelength,flux,z,dz=-0.05)

# #w,f=oz.get_range(wavelength_2,flux_2,[2600,3000])

# #plt.figure()
# #plt.plot(w,f,'k-')
# #plt.show()

# #wavelength,flux=wavelength_2,flux_2

# #print(oz.redshift_effect(x=6565,z=1.2,dz=0.1))
# #2799, dz=0.1, dx=-55
# #4863, dz=0.1, dx=-96
# #6565, dz=0.1, dx=-130




# ##### Single peak fitting of a single spectra.

# #initial_fit=[70,2800,50,1.5]
# #m=[2,1]
# #ks = [40,40,10,10]
# #intersect_limit = [2750,2850]
# #plots=[False,True,False]
# #zoom = False
# #DL=3534.4

# #oz.autorun_fit(folder_path,filename,initial_fit,m,ks,intersect_limit,plots,zoom,DL)

# ##### Generating images of the spectra

# #tf = open(folder_path+"SpARCS_IDs.txt","r")
# #contents = tf.readlines()
# #tf.close()
# #SpARCS_list=[]
# #for line in contents:
# #    SpARCS_list.append(str(line)[:-1])
# #SpARCS_desc = ['SpARCS_','.fits']
# #
# #tf = open(folder_path+"SVA1_IDs.txt","r")
# #contents = tf.readlines()
# #tf.close()
# #SVA_list=[]
# #for line in contents:
# #    SVA_list.append(str(line)[:-1])
# #SVA_desc = ['SVA1_COADD-','.fits']
# #
# #for ID in SpARCS_list:
# #    filename = SpARCS_desc[0]+ID+SpARCS_desc[1]
# #    wavelength,flux,z=oz.get_spectra(folder_path,filename,i=0,m=[7.0,7],mask_break=True,save=True)
# #
# #for ID in SVA_list:
# #    filename = SVA_desc[0]+ID+SVA_desc[1]
# #    wavelength,flux,z=oz.get_spectra(folder_path,filename,i=0,m=[7.0,7],mask_break=True,save=True)

# #for ID in ID_list:
# #    wavelength,flux=oz.get_spectra(folder_path,descriptors[0]+ID+descriptors[1],0,plot=False,save=True)

# ##### Searching the IRSA SWIRE Database for matching entries.

# #oz.swire_fields_info()
# #oz.irsa_search(folder_path,"RDZs.xlsx","RDZs_test.xlsx","elaiss1_cat_f05",1.0,[2,319],['C','D'],['F','G','H','I','J'])
# #oz.plot_on_sky(folder_path,"RDZs.xlsx",[2,319],['C','D'])

# ##### Exploring the headers of the files.

# #ID_range=[11,21]
# #descriptors = ['SVA1_COADD-','.fits']
# #cols = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O']
# #flags = ["FITSFILE","INDEX","SIMPLE","XTENSION","EXTNAME","SOURCE","RA","DEC","NDATE","Z","TOTALEXP","OBJECT","UTDATE","UTSTART","UTEND"]
# #oz.make_ID_list(folder_path, "test.txt", descriptors[0], ID_range)
# #oz.extract_reduced_headers(folder_path,"SVA1_IDs.txt","test.xlsx",descriptors,cols,flags)

# ##### Original contents of the ozdes_spectra.py file, modified by myself.

# #folder_path='C:/Users/Val/Desktop/SpARCS_2018Dec19T113006/'
# #file_name='SVA1_COADD-2937428417.fits'
# #hdr_list=fits.open(folder_path+file_name)
# #hdr=hdr_list[0].header
# #data=hdr_list[1].data
# #flux=np.array(data)
# #lam1=hdr['CRVAL1']
# #step=hdr['CDELT1']
# #z=hdr['Z']
# #print(4959.*(z+1))
# #lam=np.arange(lam1,lam1+5000.*step,step)
# #plt.plot(lam,flux)
# #plt.yscale('log')
# #plt.show()
