#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 10:33:32 2018

Add materials to smesh file.

@author: Fabricio Fernandez <fabricio.hmf@gmail.com>

"""

import sys
import os
import pandas as pd

#fn_smesh = sys.argv[1]
#fn_mat = sys.argv[2]

#--
f_in  = open('ElementSelection#1.smesh','r')
f_out = open('smesh_with_materials.smesh','w')

#--
df = pd.read_csv('materiales.csv', sep=',', header=0)
imat = df.MAT

#--
print('reading/writing points...')
label = ['SMESH.ELEMENTS\n']

while True:
    
    line = f_in.readline()
    
    if line != label[0]:
        f_out.write(line)

    else:
        f_out.write('SMESH.MATERIAL.FLAG\n\n')
        f_out.write(line)
        line = f_in.readline()
        f_out.write(line)
        break

#--
# f_out.write('SMESH.MATERIAL.FLAG\n\n')
#f_out.write('SMESH.ELEMENTS\n')
# f_in.readline()
# f_out.write(line)
#--
print('writing elements...')
count = 0
while True:

    if count==len(imat):
        break

    line = f_in.readline()
    splited = line.split()

    for i in range(len(splited)):
        
        if i!=(len(splited)-1):
            f_out.write(splited[i]+' ')
        else:
            f_out.write(' '+splited[i]+' '+str(imat[count])+'\n')
        
    #     if istr==line.split()[-1]:
    #         f_out.write(' '+str(imat[count])+'\n')

    # for istr in line.split():
        
    #     f_out.write(istr+' ')
        
    #     if istr==line.split()[-1]:
    #         f_out.write(' '+str(imat[count])+'\n')

    count = count + 1

    if line=='':
        break
#--
f_in.close()
f_out.close()