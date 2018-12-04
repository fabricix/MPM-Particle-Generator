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


fn_in = 'materiales.dat'
fn_out = 'materiales.csv'


f_in   = open(fn_in,'r')
nlines = len(f_in.readlines())
f_in.close()


#files
f_in   = open(fn_in,'r')
f_out  = open(fn_out,'w')

# read material file
for il in range(nlines):
	
	line = f_in.readline()

	for istr in line.split():
		
		f_out.write(istr)
		
		if istr !=line.split()[-1]:
			f_out.write(",")
		
		if istr==line.split()[-1]:
			f_out.write("\n")

f_in.close()
f_out.close()

# # create global coordinates
# for i in range(len(x_ds)):
	
#     # local material coords 
#     xl = x_ds[i];
#     yl = 0.0;
#     zl = z_ds[i];

#     # global coordinates of material
#     xgl = x_o[isec] + (xl*cos_phi[isec] - yl*sin_phi[isec])
#     ygl = y_o[isec] + (xl*sin_phi[isec] + yl*cos_phi[isec])
#     zgl = zl

#     # append values in global array
#     x_global.append(float(xgl))
#     y_global.append(float(ygl))
#     z_global.append(float(zgl))
#     mat_global.append(float(matid[matstr]))



# # auxiliar labels
# label = ['SMESH.ELEMENTS\n']


# print('reading ' + fn_mpm + ' file...' )

# while True:
	
#     line = f_mpm.readline()
	
#     if line!='' and line != label[0]:
#         f_out.write(line)
#     else:
#         break

# print('creating ' + fn_out + ' file...' )

# f_mpm.close()

# ## write the particles in mpm.dat

# flat_par = False
# lines = f_par.readlines()
# count = -2
# for il in lines:
	
#     line = il
#     # print(line)

#     if line == label[0]:
#         flat_par = True

#     if flat_par == True:
#         f_out.write(line)
#         count = count+1

# print('adding particles,', count,' in ' + fn_out + ' file ...' )


# f_par.close()
# f_out.close()

# ##############################################################################
# # remove auxiliar files
# # 
# ##############################################################################

# try:
#     print('removing ' + fn_mpm + ' file...')
#     os.remove(fn_mpm)
#     print('renaming  ' + fn_out + ' file to ' + fn_mpm )
#     os.rename(fn_out, fn_mpm)
# except:
#     print(' **** Not was possible to renaming and removig the files...')
 
# print('done!')