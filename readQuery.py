#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 15 15:21:37 2020

@author: val
"""

import numpy as np
import pandas as pd
import os.path

fileName = "Source_Locations.txt"
names = np.loadtxt(fileName, dtype={'names':('ID', 'RA', 'DEC'), 'formats': (np.int, np.float, np.float)})

dataLoc = "query/DESDR_"
procLoc = "processed_query.txt"


file = open(procLoc,"a+") 
 
file.write("sID\tsRA\tsDEC\tfID\tTILE\tfRA\tfDEC\tDIST\tMAGGV\tMAGGE\tMAGRV\tMAGRE\tMAGIV\tMAGIE\tBADQUAL\n")

for i in range(len(names)):
    ra = names['RA'][i]
    dec = names['DEC'][i]
    sID = names['ID'][i]
    
    print("Processing: " + str(sID))
    
    if not os.path.isfile(dataLoc + str(sID) + ".tab"):
        print("No associated entry.")
        entry=str(sID)+"\t"+str(ra)+"\t"+str(dec)+"\n"
        file.write(entry)
        continue
    
    data = pd.read_table(dataLoc + str(sID) + ".tab", delim_whitespace=True)
    
    for j in range(len(data)):
        
        distance=np.sqrt((ra-data['RA'][j])**2+(dec-data['DEC'][j])**2)
        
        entry=str(sID)+"\t"+str(ra)+"\t"+str(dec)+"\t"
        entry+=str(data['COADD_OBJECT_ID'][j])+"\t"+str(data['TILENAME'][j])+"\t"
        entry+=str(data['RA'][j])+"\t"+str(data['DEC'][j])+"\t"+str(distance)+"\t"
        entry+=str(data['MAG_AUTO_G'][j])+"\t"+str(data['MAGERR_AUTO_G'][j])+"\t"
        entry+=str(data['MAG_AUTO_R'][j])+"\t"+str(data['MAGERR_AUTO_R'][j])+"\t"
        entry+=str(data['MAG_AUTO_I'][j])+"\t"+str(data['MAGERR_AUTO_I'][j])+"\t"
    
        flags=""
        if data['FLAGS_G'][j]>3:
            flags+="G "
        if data['FLAGS_R'][j]>3:
            flags+="R "
        if data['FLAGS_I'][j]>3:
            flags+="I "
        entry+=flags+"\n"
    
        file.write(entry)
        
file.close() 
    
    
    
    
    