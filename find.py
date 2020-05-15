# -*- coding: utf-8 -*-
"""
Created on Thu May  7 16:56:18 2020

@author: ambar
"""

import sp
from astropy.io import ascii
import time
import ccdproc as ccdp
import numpy as np

a='fasfas'
print(type(a))

science=ccdp.ImageFileCollection('Corrected_Science')

#ascii.write(sc,'sc_summary.csv',format='csv',overwrite=True)

for i,j in science.hdus(return_fname=True):
    print(i)
#found=np.zeros(len(hdulist),dtype=type(a))
'''
start_time=time.time()

for imsc, imsd in enumerate(science.data()):
    f=sp.find(imsd,hmin=1000,fwhm=12)
    found[imsc]=len(f[0])

ascii.write(found,'found.csv',format='csv',overwrite=True)
'''