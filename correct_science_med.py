import glob
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.nddata import CCDData
import astropy.units as u
import ccdproc as ccdp
import os
import pathlib
from ccdproc import ImageFileCollection
from astropy.visualization import hist
import itertools
from astropy.stats import sigma_clip, mad_std
import time
import sys


print('assuming folder \'Master_Files/Windowed\' has been populated')


if not os.path.exists('Final_Science'):
    os.makedirs('Final_Science')
    print('folder \'Final_Science\' created')
if not os.path.exists('Final_Science/median'):
    os.makedirs('Final_Science/median')
    print('folder \'Final_Science/median\' created')


#read science images and find RTDATSEC
sciencecollection=ImageFileCollection('HD115709/SII',ext=1)
for sci in sciencecollection.ccds(ccd_kwargs={'unit': 'adu'}):
    sciwin=sci.meta['RTDATSEC']
    break

cmpath='Master_Files/Windowed/'
cspath='Final_Science/median/'


tmastercollection=ImageFileCollection(cmpath)


'''check shape'''
#check shape
# print(sci.data.shape)
# for master in tmastercollection.ccds():
#     print(master.data.shape)


for mbias, mbiasn in tmastercollection.ccds(imtype='mbias', combined='median',return_fname=True):
    print('using',mbiasn)
    for mflat, mflatn in tmastercollection.ccds(imtype='mflat', flatcom='median', return_fname=True):
        print('using', mflatn)
        normmed=mflat.header['normmed']
        start=time.time()
        for sci, scin in sciencecollection.ccds(return_fname=True, ccd_kwargs={'unit': 'adu'}):
            print('correcting',scin)
            subsci=ccdp.subtract_bias(sci,mbias,add_keyword={'subsci':'sigmbias'})
            corsci=ccdp.flat_correct(subsci,mflat,norm_value=float(normmed))
            corsci.write(cspath+scin,overwrite=True)

print (time.time()-start)