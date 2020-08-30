#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@created on: 19:09:20 23:06:46
@last saved: 20:02:27 23:13:23
@author: val 
"""

filepath="/Users/val/Desktop/university/"
filename="observation_tables.txt"

print("Loading "+filename+" from "+filepath)
with open(filepath+filename) as f:
    lines = f.readlines()

print("Scanning...")
targets,ras,decs=[],[],[]
for line in lines:
    entries=line.split(" ")
    if entries[0] == 'Source:':
        targets.append(entries[-1].strip("\n"))
    if entries[0] == 'RA:':
        ras.append(float(entries[-1].strip("\n")))
    if entries[0] == 'DEC:':
        decs.append(float(entries[-1].strip("\n")))
print("Entries found: "+str(len(targets)))
print("Reformatting to new txt file...")

filename="obs_table.txt"
f = open(filepath+filename, "w+") 
f.write("Target\tRA\tDEC\n")
for i in range(len(targets)):
    entry=targets[i]+"\t"
    entry+="{:.6f}".format(ras[i])+"\t"
    entry+="{:.6f}".format(decs[i])+"\n"
    f.write(entry)
f.close()

print("Completed: saved under "+filename+" in working directory.")