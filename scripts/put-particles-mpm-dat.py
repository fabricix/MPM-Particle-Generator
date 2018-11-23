#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 21:47:02 2017

Elimitate node end Elements data from mpm.dat file

@author: Fabricio Fernandez <fabricio.hmf@gmail.com>

"""

import sys
import os

#fn_mesh = sys.argv[1]
#fn_write = sys.argv[2]

#names
fn_mpm    = 'mpm.dat'
fn_out    = 'auxmpm.dat'
fn_par    = 'mpm.part'

#files
f_mpm   = open(fn_mpm,'r')
f_par   = open(fn_par,'r')
f_out   = open(fn_out,'w')

# auxiliar labels
label = ['%PARTICLES\n']

#############################################################################
# find the nodes in mesh
# create a nodecoord file
#############################################################################


## eliminate the particles in mpm.dat

print('reading ' + fn_mpm + ' file...' )

while True:
    
    line = f_mpm.readline()
    
    if line!='' and line != label[0]:
        f_out.write(line)
    else:
        break

print('creating ' + fn_out + ' file...' )

f_mpm.close()

## write the particles in mpm.dat

flat_par = False
lines = f_par.readlines()
count = -2
for il in lines:
    
    line = il
    # print(line)

    if line == label[0]:
        flat_par = True

    if flat_par == True:
        f_out.write(line)
        count = count+1

print('adding particles,', count,' in ' + fn_out + ' file ...' )


f_par.close()
f_out.close()

##############################################################################
# remove auxiliar files
# 
##############################################################################

try:
    print('removing ' + fn_mpm + ' file...')
    os.remove(fn_mpm)
    print('renaming  ' + fn_out + ' file to ' + fn_mpm )
    os.rename(fn_out, fn_mpm)
except:
    print(' **** Not was possible to renaming and removig the files...')
 
print('done!')